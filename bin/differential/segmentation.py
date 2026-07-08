from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.special import logsumexp

from .posterior import GridPosterior, normal_grid_log_mass, normalize_log_mass


@dataclass(frozen=True, slots=True)
class LengthPrior:
    """Explicit length weights with an optional inferred geometric tail."""

    log_mass: np.ndarray
    log_tail_ratio: np.ndarray | float = -np.inf

    def __post_init__(self) -> None:
        value = np.asarray(self.log_mass, dtype=float)
        ratio = np.asarray(self.log_tail_ratio, dtype=float)
        if value.ndim == 1:
            value = value[None]
        if ratio.ndim == 0:
            ratio = np.full(value.shape[0], ratio)
        weight = np.exp(value - np.max(value, axis=1, keepdims=True))
        r = np.exp(ratio)
        if np.any(r >= 1):
            raise ValueError("tail ratios must be below one")
        total = weight.sum(axis=1) + weight[:, -1] * r / (1 - r)
        with np.errstate(divide="ignore"):
            value = np.log(weight / total[:, None])
        object.__setattr__(self, "log_mass", np.ascontiguousarray(value))
        object.__setattr__(self, "log_tail_ratio", np.ascontiguousarray(ratio))

    @classmethod
    def from_log_mass(cls, log_mass: np.ndarray, infer_tail: bool = False) -> "LengthPrior":
        value = np.asarray(log_mass, dtype=float)
        rows = value[None] if value.ndim == 1 else value
        ratio = np.full(rows.shape[0], -np.inf)
        if infer_tail:
            ratio = rows[:, -1] - rows[:, -2]
            if np.any(ratio >= 0):
                raise ValueError("the inferred tail ratio must be below one")
        return cls(value, ratio)

    def for_states(self, n_state: int) -> tuple[np.ndarray, np.ndarray]:
        mass = self.log_mass
        ratio = self.log_tail_ratio
        if mass.shape[0] == 1:
            mass = np.broadcast_to(mass, (n_state, mass.shape[1]))
            ratio = np.broadcast_to(ratio, n_state)
        if mass.shape[0] != n_state:
            raise ValueError("length-prior state dimension does not match")
        return np.ascontiguousarray(mass), np.ascontiguousarray(ratio)


@dataclass(frozen=True, slots=True)
class Segmentation:
    names: tuple[str, ...]
    posterior: GridPosterior
    boundary: np.ndarray
    log_partition: np.ndarray

    def to_npz(self, path) -> None:
        np.savez_compressed(
            path, names=self.names, state_x=self.posterior.x,
            log_posterior=self.posterior.log_mass,
            boundary=self.boundary, log_partition=self.log_partition,
        )

    @classmethod
    def from_npz(cls, path) -> "Segmentation":
        with np.load(path, allow_pickle=False) as x:
            return cls(
                tuple(x["names"].tolist()),
                GridPosterior(x["state_x"], x["log_posterior"]),
                x["boundary"], x["log_partition"],
            )


def segment(
    loglik: np.ndarray,
    state_x: np.ndarray,
    names: tuple[str, ...],
    length_prior: LengthPrior,
    log_state_prior: np.ndarray | None = None,
    transition_sd: float | None = None,
    forbid_same_state: bool = True,
) -> Segmentation:
    loglik = np.asarray(loglik, dtype=float)
    state_x = np.asarray(state_x, dtype=float)
    if loglik.ndim == 2:
        loglik = loglik[None]
    log_state_prior = normalize_log_mass(log_state_prior, len(state_x))
    log_length, log_tail = length_prior.for_states(len(state_x))
    log_transition, log_stationary = _transition(
        state_x, log_state_prior, transition_sd, forbid_same_state
    )
    mean_length, log_survival, log_equilibrium, log_q = _length_terms(
        log_length, log_tail, loglik.shape[1]
    )
    stationary_mean = float(np.exp(log_stationary) @ mean_length)

    g, n, m = loglik.shape
    log_partition = np.empty(g)
    log_posterior = np.empty_like(loglik)
    boundary = np.empty((g, max(n - 1, 0)))
    for i in range(g):
        log_partition[i], log_posterior[i], boundary[i] = _partition(
            np.ascontiguousarray(loglik[i]), log_q, log_survival,
            log_equilibrium, log_transition, log_stationary, stationary_mean,
        )
    return Segmentation(names, GridPosterior(state_x, log_posterior), boundary, log_partition)


def fit_group_means(
    differential,
    length_prior: LengthPrior,
    prior_sd_floor: float | None = None,
    transition_sd: float | None = None,
    forbid_same_state: bool = False,
) -> Segmentation:
    mean = float(np.mean(differential.common_mu))
    sd = float(np.std(differential.common_mu))
    floor = prior_sd_floor or np.median(np.diff(differential.mu_x))
    prior = normal_grid_log_mass(differential.mu_x, mean, max(sd, floor))
    return segment(
        differential.loglik_mu, differential.mu_x, differential.group_names,
        length_prior, prior, transition_sd, forbid_same_state,
    )


def _transition(x, log_prior, sd, forbid_same):
    m = len(x)
    if m == 1:
        return np.zeros((1, 1)), np.zeros(1)
    log_t = np.broadcast_to(log_prior, (m, m)).copy()
    if sd is not None:
        log_t -= (x[:, None] - x[None]) ** 2 / (2 * sd * sd)
    if forbid_same:
        np.fill_diagonal(log_t, -np.inf)
    log_t -= logsumexp(log_t, axis=1, keepdims=True)
    t = np.exp(log_t)
    stationary = np.full(m, 1 / m)
    for _ in range(10000):
        new = stationary @ t
        if np.max(np.abs(new - stationary)) < 1e-14:
            break
        stationary = new
    return log_t, np.log(stationary / stationary.sum())


def _length_terms(log_mass, log_tail, n):
    m, l0 = log_mass.shape
    q = np.full((m, n), -np.inf)
    survival = np.full((m, n + 1), -np.inf)
    equilibrium = np.full((m, n + 1), -np.inf)
    mean = np.empty(m)
    for state in range(m):
        explicit = np.exp(log_mass[state])
        r = np.exp(log_tail[state])
        for length in range(1, n + 1):
            if length <= l0:
                q[state, length - 1] = log_mass[state, length - 1]
            elif r > 0:
                q[state, length - 1] = log_mass[state, -1] + (length - l0) * log_tail[state]
        for visible in range(1, n + 1):
            if visible <= l0:
                s = explicit[visible - 1:].sum() + explicit[-1] * r / (1 - r) if r < 1 else np.inf
                e = np.sum((np.arange(visible, l0 + 1) - visible + 1) * explicit[visible - 1:])
                if r > 0:
                    e += explicit[-1] * ((l0 - visible + 1) * r / (1 - r) + r / (1 - r) ** 2)
            else:
                s = explicit[-1] * r ** (visible - l0) / (1 - r) if r > 0 else 0.0
                e = explicit[-1] * r ** (visible - l0) / (1 - r) ** 2 if r > 0 else 0.0
            survival[state, visible] = np.log(s) if s > 0 else -np.inf
            equilibrium[state, visible] = np.log(e) if e > 0 else -np.inf
        base = np.sum(np.arange(1, l0 + 1) * explicit)
        if r > 0:
            base += explicit[-1] * (l0 * r / (1 - r) + r / (1 - r) ** 2)
        mean[state] = base
    return mean, survival, equilibrium, q


@njit(cache=True)
def _logadd(a, b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    x = max(a, b)
    return x + np.log(np.exp(a - x) + np.exp(b - x))


@njit(cache=True)
def _partition(emission, log_q, log_survival, log_equilibrium, log_t, log_nu, mean_length):
    n, m = emission.shape
    prefix = np.zeros((n + 1, m))
    for p in range(n):
        for s in range(m):
            prefix[p + 1, s] = prefix[p, s] + emission[p, s]

    forward = np.full((n + 1, m), -np.inf)
    incoming = np.full((n + 1, m), -np.inf)
    for end in range(1, n + 1):
        for state in range(m):
            value = log_nu[state] + log_survival[state, end] - np.log(mean_length) + prefix[end, state]
            for start in range(1, end):
                length = end - start
                term = incoming[start, state] + log_q[state, length - 1] + prefix[end, state] - prefix[start, state]
                value = _logadd(value, term)
            forward[end, state] = value
        for target in range(m):
            for source in range(m):
                incoming[end, target] = _logadd(incoming[end, target], forward[end, source] + log_t[source, target])

    reverse = np.full((n + 1, m), -np.inf)
    outgoing = np.full((n + 1, m), -np.inf)
    for start in range(n - 1, -1, -1):
        for state in range(m):
            remaining = n - start
            value = log_survival[state, remaining] + prefix[n, state] - prefix[start, state]
            for end in range(start + 1, n):
                length = end - start
                term = log_q[state, length - 1] + prefix[end, state] - prefix[start, state] + outgoing[end, state]
                value = _logadd(value, term)
            reverse[start, state] = value
        for previous in range(m):
            for target in range(m):
                outgoing[start, previous] = _logadd(outgoing[start, previous], log_t[previous, target] + reverse[start, target])

    log_z = -np.inf
    for state in range(m):
        log_z = _logadd(log_z, log_nu[state] + log_equilibrium[state, n] - np.log(mean_length) + prefix[n, state])
    for start in range(1, n):
        for state in range(m):
            log_z = _logadd(log_z, incoming[start, state] + log_survival[state, n - start] + prefix[n, state] - prefix[start, state])

    scale = np.full(m, -np.inf)
    for state in range(m):
        scale[state] = log_nu[state] + log_equilibrium[state, n] - np.log(mean_length) + prefix[n, state] - log_z
        for end in range(1, n):
            scale[state] = max(scale[state], log_nu[state] + log_survival[state, end] - np.log(mean_length) + prefix[end, state] + outgoing[end, state] - log_z)
        for start in range(1, n):
            scale[state] = max(scale[state], incoming[start, state] + log_survival[state, n - start] + prefix[n, state] - prefix[start, state] - log_z)
        for start in range(1, n - 1):
            for end in range(start + 1, n):
                scale[state] = max(scale[state], incoming[start, state] + log_q[state, end - start - 1] + prefix[end, state] - prefix[start, state] + outgoing[end, state] - log_z)

    difference = np.zeros((n + 1, m))
    for state in range(m):
        weight = np.exp(log_nu[state] + log_equilibrium[state, n] - np.log(mean_length) + prefix[n, state] - log_z - scale[state])
        difference[0, state] += weight; difference[n, state] -= weight
        for end in range(1, n):
            weight = np.exp(log_nu[state] + log_survival[state, end] - np.log(mean_length) + prefix[end, state] + outgoing[end, state] - log_z - scale[state])
            difference[0, state] += weight; difference[end, state] -= weight
        for start in range(1, n):
            weight = np.exp(incoming[start, state] + log_survival[state, n - start] + prefix[n, state] - prefix[start, state] - log_z - scale[state])
            difference[start, state] += weight; difference[n, state] -= weight
        for start in range(1, n - 1):
            for end in range(start + 1, n):
                value = incoming[start, state] + log_q[state, end - start - 1] + prefix[end, state] - prefix[start, state] + outgoing[end, state] - log_z
                weight = np.exp(value - scale[state])
                difference[start, state] += weight; difference[end, state] -= weight

    log_posterior = np.full((n, m), -np.inf)
    running = np.zeros(m)
    for p in range(n):
        normalizer = -np.inf
        for state in range(m):
            running[state] += difference[p, state]
            if running[state] > 0:
                log_posterior[p, state] = np.log(running[state]) + scale[state]
            normalizer = _logadd(normalizer, log_posterior[p, state])
        for state in range(m):
            log_posterior[p, state] -= normalizer

    boundary = np.empty(max(n - 1, 0))
    for b in range(1, n):
        value = -np.inf
        for source in range(m):
            for target in range(m):
                value = _logadd(value, forward[b, source] + log_t[source, target] + reverse[b, target])
        boundary[b - 1] = min(1.0, np.exp(value - log_z))
    return log_z, log_posterior, boundary
