from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from .config import DEFAULT_ETA_SEGMENTATION, EtaSegmentationConfig
from .posterior import GridPosterior, normalize_log_mass, spike_slab_log_mass
from .segmentation import LengthPrior, segment
from .variance_ratio import VarianceRatioLikelihood


@dataclass(frozen=True, slots=True)
class EtaSegmentation:
    group_names: tuple[str, ...]
    eta_x: np.ndarray
    mu0: GridPosterior
    icc: GridPosterior
    boundary: np.ndarray
    log_partition: float
    log_mu0_prior: np.ndarray
    log_eta_prior: np.ndarray

    def to_npz(self, path) -> None:
        np.savez_compressed(
            path, group_names=self.group_names, eta_x=self.eta_x,
            mu0_x=self.mu0.x, log_mu0=self.mu0.log_mass,
            icc_x=self.icc.x, log_icc=self.icc.log_mass,
            boundary=self.boundary, log_partition=self.log_partition,
            log_mu0_prior=self.log_mu0_prior, log_eta_prior=self.log_eta_prior,
        )

    @classmethod
    def from_npz(cls, path) -> "EtaSegmentation":
        with np.load(path, allow_pickle=False) as x:
            return cls(
                tuple(x["group_names"].tolist()), x["eta_x"],
                GridPosterior(x["mu0_x"], x["log_mu0"]),
                GridPosterior(x["icc_x"], x["log_icc"]),
                x["boundary"], float(x["log_partition"]),
                x["log_mu0_prior"], x["log_eta_prior"],
            )


def eta_log_prior(eta_x: np.ndarray, config: EtaSegmentationConfig) -> np.ndarray:
    spike = np.flatnonzero(np.isneginf(eta_x))
    if spike.size:
        return spike_slab_log_mass(
            eta_x, int(spike[0]), config.consistent_mass,
            config.slab_mean, config.slab_sd,
        )
    from .posterior import normal_grid_log_mass
    return normal_grid_log_mass(eta_x, config.slab_mean, config.slab_sd)


def fit_eta_segmentation(
    likelihood: VarianceRatioLikelihood,
    length_prior: LengthPrior,
    config: EtaSegmentationConfig = DEFAULT_ETA_SEGMENTATION,
    log_mu0_prior: np.ndarray | None = None,
    log_eta_prior: np.ndarray | None = None,
) -> EtaSegmentation:
    log_mu0_prior = normalize_log_mass(log_mu0_prior, len(likelihood.mu0_x))
    log_eta_prior = eta_log_prior(likelihood.eta_x, config) if log_eta_prior is None else normalize_log_mass(log_eta_prior, len(likelihood.eta_x))
    emission = logsumexp(likelihood.loglik + log_mu0_prior[None, :, None], axis=1)
    base = segment(
        emission, likelihood.icc_x, ("eta",), length_prior,
        log_eta_prior, config.transition_sd, config.forbid_same_state,
    )
    log_icc = base.posterior.log_mass[0]
    log_mu0 = _reconstruct_mu0(
        likelihood.loglik, log_mu0_prior, emission, log_icc
    )
    return EtaSegmentation(
        likelihood.group_names, likelihood.eta_x,
        GridPosterior(likelihood.mu0_x, log_mu0),
        GridPosterior(likelihood.icc_x, log_icc),
        base.boundary[0], float(base.log_partition[0]),
        log_mu0_prior, log_eta_prior,
    )


def _reconstruct_mu0(loglik, log_mu0_prior, emission, log_state):
    conditional = loglik + log_mu0_prior[None, :, None] - emission[:, None, :]
    return logsumexp(log_state[:, None, :] + conditional, axis=2)
