from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from .config import (
    DEFAULT_COEFFICIENT,
    DEFAULT_COEFFICIENT_SEGMENTATION,
    CoefficientConfig,
    CoefficientSegmentationConfig,
)
from .differential import Differential
from .integration import (
    coefficient_target_terms,
    leave_one_out_reference,
    variance_ratio_group_terms,
)
from .posterior import GridPosterior, normalize_log_mass, spike_slab_log_mass
from .segmentation import LengthPrior, Segmentation, segment
from .variance_ratio import VarianceRatioLikelihood


@dataclass(frozen=True, slots=True)
class CoefficientLikelihood:
    group_names: tuple[str, ...]
    z_x: np.ndarray
    loglik: np.ndarray

    def posterior(self, log_z_prior: np.ndarray | None = None) -> GridPosterior:
        prior = normalize_log_mass(log_z_prior, len(self.z_x))
        return GridPosterior(self.z_x, self.loglik + prior[None, None])

    def to_npz(self, path) -> None:
        np.savez_compressed(path, group_names=self.group_names, z_x=self.z_x, loglik=self.loglik)

    @classmethod
    def from_npz(cls, path) -> "CoefficientLikelihood":
        with np.load(path, allow_pickle=False) as x:
            return cls(tuple(x["group_names"].tolist()), x["z_x"], x["loglik"])


class CoefficientModel:
    def __init__(self, config: CoefficientConfig = DEFAULT_COEFFICIENT):
        self.config = config

    def fit(
        self,
        differential: Differential,
        variance_ratio: VarianceRatioLikelihood,
        log_mu0_prior: np.ndarray | None = None,
        log_eta_prior: np.ndarray | None = None,
    ) -> CoefficientLikelihood:
        mu_prior = normalize_log_mass(log_mu0_prior, len(variance_ratio.mu0_x))
        eta_prior = normalize_log_mass(log_eta_prior, len(variance_ratio.eta_x))
        group_terms = variance_ratio_group_terms(
            differential, variance_ratio.mu0_x, variance_ratio.eta_x,
            self.config.method,
        )
        reference = leave_one_out_reference(group_terms, eta_prior)
        target = coefficient_target_terms(
            differential, variance_ratio.mu0_x,
            self.config.z_x(), self.config.method,
        )
        loglik = logsumexp(
            reference[:, :, :, None] + target + mu_prior[None, None, :, None],
            axis=2,
        )
        return CoefficientLikelihood(
            differential.group_names, self.config.z_x(), loglik
        )


def z_log_prior(z_x: np.ndarray, config: CoefficientSegmentationConfig) -> np.ndarray:
    zero = int(np.argmin(np.abs(z_x)))
    return spike_slab_log_mass(z_x, zero, config.zero_mass, 0.0, config.slab_sd)


def fit_coefficient_segmentation(
    likelihood: CoefficientLikelihood,
    length_prior: LengthPrior,
    config: CoefficientSegmentationConfig = DEFAULT_COEFFICIENT_SEGMENTATION,
    log_z_prior: np.ndarray | None = None,
) -> Segmentation:
    prior = z_log_prior(likelihood.z_x, config) if log_z_prior is None else normalize_log_mass(log_z_prior, len(likelihood.z_x))
    return segment(
        likelihood.loglik, likelihood.z_x, likelihood.group_names,
        length_prior, prior, config.transition_sd, config.forbid_same_state,
    )
