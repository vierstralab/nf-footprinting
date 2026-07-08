from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from .config import DEFAULT_VARIANCE_RATIO, VarianceRatioConfig
from .differential import Differential
from .integration import variance_ratio_group_terms
from .posterior import GridPosterior, normalize_log_mass


@dataclass(frozen=True, slots=True)
class VarianceRatioPosterior:
    mu0: GridPosterior
    icc: GridPosterior


@dataclass(frozen=True, slots=True)
class VarianceRatioLikelihood:
    group_names: tuple[str, ...]
    mu0_x: np.ndarray
    eta_x: np.ndarray
    loglik: np.ndarray

    @property
    def ratio_x(self) -> np.ndarray:
        return np.where(np.isneginf(self.eta_x), 0.0, np.exp(self.eta_x))

    @property
    def icc_x(self) -> np.ndarray:
        ratio = self.ratio_x
        return ratio / (1 + ratio)

    def posterior(
        self,
        log_mu0_prior: np.ndarray | None = None,
        log_eta_prior: np.ndarray | None = None,
    ) -> VarianceRatioPosterior:
        joint = (
            self.loglik
            + normalize_log_mass(log_mu0_prior, len(self.mu0_x))[None, :, None]
            + normalize_log_mass(log_eta_prior, len(self.eta_x))[None, None]
        )
        joint -= logsumexp(joint, axis=(1, 2), keepdims=True)
        return VarianceRatioPosterior(
            GridPosterior(self.mu0_x, logsumexp(joint, axis=2)),
            GridPosterior(self.icc_x, logsumexp(joint, axis=1)),
        )

    def to_npz(self, path) -> None:
        np.savez_compressed(
            path, group_names=self.group_names, mu0_x=self.mu0_x,
            eta_x=self.eta_x, loglik=self.loglik,
        )

    @classmethod
    def from_npz(cls, path) -> "VarianceRatioLikelihood":
        with np.load(path, allow_pickle=False) as x:
            return cls(tuple(x["group_names"].tolist()), x["mu0_x"], x["eta_x"], x["loglik"])


class VarianceRatioModel:
    def __init__(self, config: VarianceRatioConfig = DEFAULT_VARIANCE_RATIO):
        self.config = config

    def fit(self, differential: Differential) -> VarianceRatioLikelihood:
        mu0_x = differential.mu_x
        eta_x = self.config.eta_x()
        group_terms = variance_ratio_group_terms(
            differential, mu0_x, eta_x,
            self.config.method, self.config.variance_floor,
        )
        return VarianceRatioLikelihood(
            differential.group_names, mu0_x, eta_x, group_terms.sum(axis=0)
        )
