from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp, ndtr


@dataclass(frozen=True, slots=True)
class GridPosterior:
    """Normalized posterior on a one-dimensional grid; state is the last axis."""

    x: np.ndarray
    log_mass: np.ndarray

    def __post_init__(self) -> None:
        x = np.asarray(self.x, dtype=float)
        log_mass = np.asarray(self.log_mass, dtype=float)
        if log_mass.shape[-1] != x.size:
            raise ValueError("posterior state axis does not match the grid")
        log_mass = log_mass - logsumexp(log_mass, axis=-1, keepdims=True)
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "log_mass", np.ascontiguousarray(log_mass))

    @classmethod
    def from_log_weight(cls, x: np.ndarray, log_weight: np.ndarray) -> "GridPosterior":
        return cls(x, log_weight)

    @property
    def mass(self) -> np.ndarray:
        return np.exp(self.log_mass)

    @property
    def mean(self) -> np.ndarray:
        return self.mass @ self.x

    @property
    def map(self) -> np.ndarray:
        return self.x[np.argmax(self.log_mass, axis=-1)]

    def quantile(self, probability: float) -> np.ndarray:
        if not 0 <= probability <= 1:
            raise ValueError("probability must lie in [0, 1]")
        return self.x[np.argmax(np.cumsum(self.mass, axis=-1) >= probability, axis=-1)]

    def prob_gt(self, threshold: float) -> np.ndarray:
        return self._prob(self.x > threshold)

    def prob_ge(self, threshold: float) -> np.ndarray:
        return self._prob(self.x >= threshold)

    def prob_lt(self, threshold: float) -> np.ndarray:
        return self._prob(self.x < threshold)

    def prob_le(self, threshold: float) -> np.ndarray:
        return self._prob(self.x <= threshold)

    def prob_abs_gt(self, threshold: float) -> np.ndarray:
        return self._prob(np.abs(self.x) > threshold)

    def prob_abs_lt(self, threshold: float) -> np.ndarray:
        return self._prob(np.abs(self.x) < threshold)

    def with_grid(self, x: np.ndarray) -> "GridPosterior":
        return GridPosterior(x, self.log_mass)

    def _prob(self, mask: np.ndarray) -> np.ndarray:
        if not np.any(mask):
            return np.zeros(self.log_mass.shape[:-1])
        return np.exp(logsumexp(self.log_mass[..., mask], axis=-1))


def normalize_log_mass(value: np.ndarray | None, size: int) -> np.ndarray:
    if value is None:
        return np.full(size, -np.log(size))
    value = np.asarray(value, dtype=float)
    return value - logsumexp(value)


def normal_grid_log_mass(x: np.ndarray, mean: float, sd: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == 1:
        return np.zeros(1)
    bounds = np.r_[-np.inf, (x[:-1] + x[1:]) / 2, np.inf]
    mass = np.maximum(np.diff(ndtr((bounds - mean) / sd)), np.finfo(float).tiny)
    return np.log(mass) - np.log(mass.sum())


def spike_slab_log_mass(
    x: np.ndarray,
    spike_index: int,
    spike_mass: float,
    slab_mean: float,
    slab_sd: float,
) -> np.ndarray:
    finite = np.ones(len(x), dtype=bool)
    finite[spike_index] = False
    slab = np.exp(normal_grid_log_mass(np.asarray(x)[finite], slab_mean, slab_sd))
    mass = np.zeros(len(x))
    mass[spike_index] = spike_mass
    mass[finite] = (1 - spike_mass) * slab
    with np.errstate(divide="ignore"):
        return np.log(mass)
