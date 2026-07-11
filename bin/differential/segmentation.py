from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.special import logsumexp

from .posterior import GridPosterior


_INDEPENDENT = 0
_INDEPENDENT_NO_SELF = 1
_SPARSE = 2


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
    log_state_prior: np.ndarray,
    transition_sd: float | None = None,
    forbid_same_state: bool = True,
) -> Segmentation:
    """Segment one or more tracks on a shared state grid.

    ``log_state_prior`` may have shape ``state`` or ``track x state``.
    Independent transitions are evaluated in linear time in the number of
    states. Gaussian attraction transitions retain only numerically nonzero
    entries, so narrow kernels are evaluated as banded sparse transitions.
    """
    loglik = np.asarray(loglik, dtype=float)
    state_x = np.asarray(state_x, dtype=float)
    if loglik.ndim == 2:
        loglik = loglik[None]
    if loglik.ndim != 3:
        raise ValueError("loglik must have shape position x state or track x position x state")

    g, n, m = loglik.shape
    if state_x.shape != (m,):
        raise ValueError("state grid does not match the likelihood state axis")
    if len(names) != g:
        raise ValueError("names must match the number of tracks")

    log_state_prior = _normalize_track_prior(log_state_prior, g, m)
    log_length, log_tail = length_prior.for_states(m)
    mean_length, log_survival, log_equilibrium = _length_terms(
        log_length, log_tail, n
    )

    if transition_sd is None:
        mode = _INDEPENDENT_NO_SELF if forbid_same_state else _INDEPENDENT
        log_nu = np.empty_like(log_state_prior)
        if mode == _INDEPENDENT:
            log_nu[:] = log_state_prior
        else:
            p = np.exp(log_state_prior)
            if np.any(p >= 1):
                raise ValueError("forbid_same_state requires at least two states with positive prior mass")
            log_nu[:] = log_state_prior + np.log1p(-p)
            log_nu -= logsumexp(log_nu, axis=1, keepdims=True)
        stationary_mean = np.exp(log_nu) @ mean_length
        log_partition = np.empty(g)
        log_posterior = np.empty_like(loglik)
        boundary = np.empty((g, max(n - 1, 0)))
        empty_i = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0)

        def run_track(i):
            return _partition(
                np.ascontiguousarray(loglik[i]), log_length, log_tail,
                log_survival, log_equilibrium, log_state_prior[i], mode,
                empty_i, empty_i, empty_f, empty_i, empty_i, empty_f,
                log_nu[i], stationary_mean[i],
            )

        results = [run_track(0)]
        if g > 1:
            with ThreadPoolExecutor(max_workers=min(g - 1, 32)) as pool:
                results.extend(pool.map(run_track, range(1, g)))
        for i, (z, posterior, borders) in enumerate(results):
            log_partition[i] = z
            log_posterior[i] = posterior
            boundary[i] = borders
    else:
        log_partition = np.empty(g)
        log_posterior = np.empty_like(loglik)
        boundary = np.empty((g, max(n - 1, 0)))
        for i in range(g):
            transition = _sparse_transition(
                state_x, log_state_prior[i], transition_sd, forbid_same_state
            )
            stationary_mean = float(np.exp(transition[-1]) @ mean_length)
            log_partition[i], log_posterior[i], boundary[i] = _partition(
                np.ascontiguousarray(loglik[i]), log_length, log_tail,
                log_survival, log_equilibrium, log_state_prior[i],
                _SPARSE, *transition[:-1], transition[-1], stationary_mean,
            )

    return Segmentation(
        names, GridPosterior(state_x, log_posterior), boundary, log_partition
    )


def _normalize_track_prior(value: np.ndarray, n_track: int, n_state: int) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    if value.ndim == 1:
        if value.shape != (n_state,):
            raise ValueError(f"prior must have shape ({n_state},) or ({n_track}, {n_state})")
        value = np.broadcast_to(value, (n_track, n_state)).copy()
    elif value.shape != (n_track, n_state):
        raise ValueError(f"prior must have shape ({n_state},) or ({n_track}, {n_state})")
    value -= logsumexp(value, axis=1, keepdims=True)
    return np.ascontiguousarray(value)


def _sparse_transition(x, log_prior, sd, forbid_same):
    """Return row/column sparse forms and the stationary state mass."""
    m = len(x)
    if m == 1:
        empty_i = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=float)
        ptr = np.zeros(2, dtype=np.int64)
        return ptr, empty_i, empty_f, ptr, empty_i, empty_f, np.zeros(1)

    log_t = np.broadcast_to(log_prior, (m, m)).copy()
    log_t -= (x[:, None] - x[None]) ** 2 / (2 * sd * sd)
    if forbid_same:
        np.fill_diagonal(log_t, -np.inf)
    log_t -= logsumexp(log_t, axis=1, keepdims=True)

    # Conversion to probabilities removes only transitions below float64's
    # representable range and gives a compact band for narrow attraction.
    t = np.exp(log_t)
    t /= t.sum(axis=1, keepdims=True)
    stationary = np.full(m, 1 / m)
    for _ in range(10000):
        new = stationary @ t
        if np.max(np.abs(new - stationary)) < 1e-14:
            stationary = new
            break
        stationary = new
    log_nu = np.log(stationary / stationary.sum())

    rows, cols = np.nonzero(t)
    values = np.log(t[rows, cols])
    row_ptr = np.zeros(m + 1, dtype=np.int64)
    np.add.at(row_ptr, rows + 1, 1)
    np.cumsum(row_ptr, out=row_ptr)

    order = np.lexsort((rows, cols))
    c_rows = rows[order].astype(np.int64)
    c_cols = cols[order]
    c_values = values[order]
    col_ptr = np.zeros(m + 1, dtype=np.int64)
    np.add.at(col_ptr, c_cols + 1, 1)
    np.cumsum(col_ptr, out=col_ptr)

    return (
        np.ascontiguousarray(row_ptr),
        np.ascontiguousarray(cols.astype(np.int64)),
        np.ascontiguousarray(values),
        np.ascontiguousarray(col_ptr),
        np.ascontiguousarray(c_rows),
        np.ascontiguousarray(c_values),
        np.ascontiguousarray(log_nu),
    )


def _length_terms(log_mass, log_tail, n):
    m, l0 = log_mass.shape
    survival = np.full((m, n + 1), -np.inf)
    equilibrium = np.full((m, n + 1), -np.inf)
    mean = np.empty(m)
    for state in range(m):
        explicit = np.exp(log_mass[state])
        r = np.exp(log_tail[state])
        for visible in range(1, n + 1):
            if visible <= l0:
                s = explicit[visible - 1:].sum() + explicit[-1] * r / (1 - r)
                e = np.sum((np.arange(visible, l0 + 1) - visible + 1) * explicit[visible - 1:])
                if r > 0:
                    e += explicit[-1] * (
                        (l0 - visible + 1) * r / (1 - r)
                        + r / (1 - r) ** 2
                    )
            else:
                s = explicit[-1] * r ** (visible - l0) / (1 - r) if r > 0 else 0.0
                e = explicit[-1] * r ** (visible - l0) / (1 - r) ** 2 if r > 0 else 0.0
            survival[state, visible] = np.log(s) if s > 0 else -np.inf
            equilibrium[state, visible] = np.log(e) if e > 0 else -np.inf
        base = np.sum(np.arange(1, l0 + 1) * explicit)
        if r > 0:
            base += explicit[-1] * (
                l0 * r / (1 - r) + r / (1 - r) ** 2
            )
        mean[state] = base
    return mean, survival, equilibrium


@njit(cache=True, inline="always")
def _logadd(a, b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    x = max(a, b)
    return x + np.log(np.exp(a - x) + np.exp(b - x))


@njit(cache=True, inline="always")
def _logsum(values):
    result = -np.inf
    for i in range(values.size):
        result = _logadd(result, values[i])
    return result


@njit(cache=True)
def _excluding_logsum(values, out, prefix, suffix):
    m = values.size
    prefix[0] = -np.inf
    for i in range(m):
        prefix[i + 1] = _logadd(prefix[i], values[i])
    suffix[m] = -np.inf
    for i in range(m - 1, -1, -1):
        suffix[i] = _logadd(values[i], suffix[i + 1])
    for i in range(m):
        out[i] = _logadd(prefix[i], suffix[i + 1])


@njit(cache=True)
def _forward_transition(values, out, mode, log_prior, log_denom,
                        col_ptr, row_idx, col_logp, scratch, prefix, suffix):
    m = values.size
    if mode == _INDEPENDENT:
        total = _logsum(values)
        for target in range(m):
            out[target] = total + log_prior[target]
    elif mode == _INDEPENDENT_NO_SELF:
        for source in range(m):
            scratch[source] = values[source] - log_denom[source]
        _excluding_logsum(scratch, out, prefix, suffix)
        for target in range(m):
            out[target] += log_prior[target]
    else:
        for target in range(m):
            value = -np.inf
            for j in range(col_ptr[target], col_ptr[target + 1]):
                value = _logadd(value, values[row_idx[j]] + col_logp[j])
            out[target] = value


@njit(cache=True)
def _reverse_transition(values, out, mode, log_prior, log_denom,
                        row_ptr, col_idx, row_logp, scratch, prefix, suffix):
    m = values.size
    if mode == _INDEPENDENT:
        for target in range(m):
            scratch[target] = log_prior[target] + values[target]
        total = _logsum(scratch)
        for source in range(m):
            out[source] = total
    elif mode == _INDEPENDENT_NO_SELF:
        for target in range(m):
            scratch[target] = log_prior[target] + values[target]
        _excluding_logsum(scratch, out, prefix, suffix)
        for source in range(m):
            out[source] -= log_denom[source]
    else:
        for source in range(m):
            value = -np.inf
            for j in range(row_ptr[source], row_ptr[source + 1]):
                value = _logadd(value, row_logp[j] + values[col_idx[j]])
            out[source] = value


@njit(cache=True)
def _boundary_transition(forward, reverse, mode, log_prior, log_denom,
                         row_ptr, col_idx, row_logp, scratch, excluded,
                         prefix, suffix):
    m = forward.size
    if mode == _INDEPENDENT:
        for target in range(m):
            scratch[target] = log_prior[target] + reverse[target]
        return _logsum(forward) + _logsum(scratch)
    if mode == _INDEPENDENT_NO_SELF:
        for target in range(m):
            scratch[target] = log_prior[target] + reverse[target]
        _excluding_logsum(scratch, excluded, prefix, suffix)
        value = -np.inf
        for source in range(m):
            value = _logadd(
                value, forward[source] - log_denom[source] + excluded[source]
            )
        return value
    value = -np.inf
    for source in range(m):
        for j in range(row_ptr[source], row_ptr[source + 1]):
            value = _logadd(
                value,
                forward[source] + row_logp[j] + reverse[col_idx[j]],
            )
    return value


@njit(cache=True, nogil=True)
def _partition(emission, log_q, log_tail, log_survival, log_equilibrium,
               log_prior, mode, row_ptr, col_idx, row_logp,
               col_ptr, row_idx, col_logp, log_nu, mean_length):
    n, m = emission.shape
    l0 = log_q.shape[1]
    log_mean = np.log(mean_length)
    log_denom = np.empty(m)
    for state in range(m):
        p = np.exp(log_prior[state])
        log_denom[state] = np.log1p(-p) if mode == _INDEPENDENT_NO_SELF else 0.0

    prefix = np.zeros((n + 1, m))
    for p in range(n):
        for s in range(m):
            prefix[p + 1, s] = prefix[p, s] + emission[p, s]

    scratch = np.empty(m)
    excluded = np.empty(m)
    work_prefix = np.empty(m + 1)
    work_suffix = np.empty(m + 1)

    forward = np.full((n + 1, m), -np.inf)
    incoming = np.full((n + 1, m), -np.inf)
    tail_acc = np.full(m, -np.inf)
    for end in range(1, n + 1):
        eligible = end - l0 - 1
        for state in range(m):
            if eligible >= 1 and log_tail[state] > -np.inf:
                term = (
                    incoming[eligible, state]
                    - prefix[eligible, state]
                    - eligible * log_tail[state]
                )
                tail_acc[state] = _logadd(tail_acc[state], term)

            value = (
                log_nu[state] + log_survival[state, end]
                - log_mean + prefix[end, state]
            )
            max_length = min(l0, end - 1)
            for length in range(1, max_length + 1):
                start = end - length
                if start >= 1 and log_q[state, length - 1] > -np.inf:
                    value = _logadd(
                        value,
                        incoming[start, state] + log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state],
                    )
            if tail_acc[state] > -np.inf:
                tail_const = log_q[state, l0 - 1] - l0 * log_tail[state]
                value = _logadd(
                    value,
                    prefix[end, state] + end * log_tail[state]
                    + tail_const + tail_acc[state],
                )
            forward[end, state] = value

        _forward_transition(
            forward[end], incoming[end], mode, log_prior, log_denom,
            col_ptr, row_idx, col_logp, scratch, work_prefix, work_suffix,
        )

    reverse = np.full((n + 1, m), -np.inf)
    outgoing = np.full((n + 1, m), -np.inf)
    tail_acc[:] = -np.inf
    for start in range(n - 1, -1, -1):
        eligible = start + l0 + 1
        for state in range(m):
            if eligible < n and log_tail[state] > -np.inf:
                term = (
                    prefix[eligible, state] + outgoing[eligible, state]
                    + eligible * log_tail[state]
                )
                tail_acc[state] = _logadd(tail_acc[state], term)

            remaining = n - start
            value = (
                log_survival[state, remaining]
                + prefix[n, state] - prefix[start, state]
            )
            max_length = min(l0, n - start - 1)
            for length in range(1, max_length + 1):
                end = start + length
                if log_q[state, length - 1] > -np.inf:
                    value = _logadd(
                        value,
                        log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state]
                        + outgoing[end, state],
                    )
            if tail_acc[state] > -np.inf:
                tail_const = log_q[state, l0 - 1] - l0 * log_tail[state]
                value = _logadd(
                    value,
                    tail_const - prefix[start, state]
                    - start * log_tail[state] + tail_acc[state],
                )
            reverse[start, state] = value

        _reverse_transition(
            reverse[start], outgoing[start], mode, log_prior, log_denom,
            row_ptr, col_idx, row_logp, scratch, work_prefix, work_suffix,
        )

    log_z = -np.inf
    for state in range(m):
        log_z = _logadd(
            log_z,
            log_nu[state] + log_equilibrium[state, n]
            - log_mean + prefix[n, state],
        )
    for start in range(1, n):
        for state in range(m):
            log_z = _logadd(
                log_z,
                incoming[start, state] + log_survival[state, n - start]
                + prefix[n, state] - prefix[start, state],
            )

    scale = np.full(m, -np.inf)
    for state in range(m):
        scale[state] = (
            log_nu[state] + log_equilibrium[state, n]
            - log_mean + prefix[n, state] - log_z
        )
        for end in range(1, n):
            scale[state] = max(
                scale[state],
                log_nu[state] + log_survival[state, end]
                - log_mean + prefix[end, state]
                + outgoing[end, state] - log_z,
            )
        for start in range(1, n):
            scale[state] = max(
                scale[state],
                incoming[start, state] + log_survival[state, n - start]
                + prefix[n, state] - prefix[start, state] - log_z,
            )
        for start in range(1, n - 1):
            max_length = min(l0, n - start - 1)
            for length in range(1, max_length + 1):
                end = start + length
                if log_q[state, length - 1] > -np.inf:
                    scale[state] = max(
                        scale[state],
                        incoming[start, state] + log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state]
                        + outgoing[end, state] - log_z,
                    )

        if log_tail[state] > -np.inf:
            best_start = -np.inf
            tail_const = log_q[state, l0 - 1] - l0 * log_tail[state]
            for end in range(1, n):
                new_start = end - l0 - 1
                if new_start >= 1:
                    best_start = max(
                        best_start,
                        incoming[new_start, state]
                        - prefix[new_start, state]
                        - new_start * log_tail[state],
                    )
                if best_start > -np.inf:
                    scale[state] = max(
                        scale[state],
                        tail_const + best_start + prefix[end, state]
                        + outgoing[end, state] + end * log_tail[state]
                        - log_z,
                    )

    difference = np.zeros((n + 1, m))
    suffix_tail = np.full(n + 1, -np.inf)
    for state in range(m):
        weight = np.exp(
            log_nu[state] + log_equilibrium[state, n]
            - log_mean + prefix[n, state] - log_z - scale[state]
        )
        difference[0, state] += weight
        difference[n, state] -= weight

        for end in range(1, n):
            weight = np.exp(
                log_nu[state] + log_survival[state, end]
                - log_mean + prefix[end, state] + outgoing[end, state]
                - log_z - scale[state]
            )
            difference[0, state] += weight
            difference[end, state] -= weight

        for start in range(1, n):
            weight = np.exp(
                incoming[start, state] + log_survival[state, n - start]
                + prefix[n, state] - prefix[start, state]
                - log_z - scale[state]
            )
            difference[start, state] += weight
            difference[n, state] -= weight

        for start in range(1, n - 1):
            max_length = min(l0, n - start - 1)
            for length in range(1, max_length + 1):
                end = start + length
                if log_q[state, length - 1] > -np.inf:
                    value = (
                        incoming[start, state] + log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state]
                        + outgoing[end, state] - log_z
                    )
                    weight = np.exp(value - scale[state])
                    difference[start, state] += weight
                    difference[end, state] -= weight

        if log_tail[state] > -np.inf:
            tail_const = log_q[state, l0 - 1] - l0 * log_tail[state] - log_z
            suffix_tail[:] = -np.inf
            for end in range(n - 1, 0, -1):
                value = (
                    prefix[end, state] + outgoing[end, state]
                    + end * log_tail[state]
                )
                suffix_tail[end] = _logadd(value, suffix_tail[end + 1])

            for start in range(1, n - l0 - 1):
                first_end = start + l0 + 1
                value = (
                    tail_const + incoming[start, state]
                    - prefix[start, state] - start * log_tail[state]
                    + suffix_tail[first_end]
                )
                difference[start, state] += np.exp(value - scale[state])

            prefix_tail = -np.inf
            for end in range(1, n):
                new_start = end - l0 - 1
                if new_start >= 1:
                    value = (
                        incoming[new_start, state]
                        - prefix[new_start, state]
                        - new_start * log_tail[state]
                    )
                    prefix_tail = _logadd(prefix_tail, value)
                if prefix_tail > -np.inf:
                    value = (
                        tail_const + prefix[end, state]
                        + outgoing[end, state] + end * log_tail[state]
                        + prefix_tail
                    )
                    difference[end, state] -= np.exp(value - scale[state])

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
        value = _boundary_transition(
            forward[b], reverse[b], mode, log_prior, log_denom,
            row_ptr, col_idx, row_logp, scratch, excluded,
            work_prefix, work_suffix,
        )
        boundary[b - 1] = min(1.0, np.exp(value - log_z))
    return log_z, log_posterior, boundary
