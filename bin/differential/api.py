from genome_tools.plotting.modular_plot.api import PlotDataLoader

from .coefficients import CoefficientModel, fit_coefficient_segmentation
from .config import (
    DEFAULT_COEFFICIENT,
    DEFAULT_COEFFICIENT_SEGMENTATION,
    DEFAULT_DIFFERENTIAL,
    DEFAULT_ETA_SEGMENTATION,
    DEFAULT_MEAN_SEGMENTATION,
    DEFAULT_VARIANCE_RATIO,
)
from .differential import DifferentialModel, fit_group_means_segmentation
from .eta import fit_eta_segmentation
from .group_counts import infer_kdev, infer_kfp
from .theta import fit_theta_likelihood, fit_theta_segmentation
from .variance_ratio import VarianceRatioModel, fit_mu0_segmentation


class DifferentialLoader(PlotDataLoader):
    def _load(self, data, selected_groups=None, config=DEFAULT_DIFFERENTIAL):
        data.differential = DifferentialModel(config).fit(
            data.groups_data,
            data.obs,
            data.exp,
            data.disp_models,
            selected_groups,
        )
        return data


class GroupMeanSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_MEAN_SEGMENTATION,
        log_mu_prior=None,
    ):
        data.segmentation = fit_group_means_segmentation(
            data.differential,
            length_prior,
            config,
            log_mu_prior,
        )
        return data


class VarianceRatioLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_VARIANCE_RATIO,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.variance_ratio = VarianceRatioModel(config).fit(
            data.differential,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class EtaSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_ETA_SEGMENTATION,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.eta_segmentation = fit_eta_segmentation(
            data.variance_ratio,
            length_prior,
            config,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class Mu0SegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_MEAN_SEGMENTATION,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.mu0_segmentation = fit_mu0_segmentation(
            data.variance_ratio,
            length_prior,
            config,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class CoefficientLikelihoodLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_COEFFICIENT,
        log_mu0_prior=None,
        log_eta_prior=None,
        log_z_prior=None,
    ):
        data.group_coefficient_likelihood = CoefficientModel(config).fit(
            data.differential,
            data.variance_ratio,
            log_mu0_prior,
            log_eta_prior,
            log_z_prior,
        )
        return data


class CoefficientSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_COEFFICIENT_SEGMENTATION,
        log_z_prior=None,
    ):
        data.group_coefficient_segmentation = fit_coefficient_segmentation(
            data.group_coefficient_likelihood,
            length_prior,
            config,
            log_z_prior,
        )
        return data


class FootprintCountLoader(PlotDataLoader):
    def _load(
        self,
        data,
        threshold=0.0,
        source="segmented",
        log_mu_prior=None,
    ):
        if source == "segmented":
            mu_posterior = data.segmentation.posterior
        elif source == "pointwise":
            mu_posterior = data.differential.posterior(log_mu_prior)
        else:
            raise ValueError("source must be 'segmented' or 'pointwise'")
        data.kfp = infer_kfp(mu_posterior, threshold)
        return data


class DeviantGroupCountLoader(PlotDataLoader):
    def _load(
        self,
        data,
        threshold=1.0,
        variance_floor=None,
        position_chunk_size=64,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.kdev = infer_kdev(
            data.differential,
            data.variance_ratio,
            threshold=threshold,
            variance_floor=variance_floor,
            position_chunk_size=position_chunk_size,
            log_mu0_prior=log_mu0_prior,
            log_eta_prior=log_eta_prior,
        )
        return data


class ThetaLoader(PlotDataLoader):
    def _load(
        self,
        data,
        mode="group_informed",
        theta_x=None,
        theta_step=0.1,
        theta_tail_sd=5.0,
        position_chunk_size=8,
        storage_dtype="float32",
        log_mu_prior=None,
    ):
        data.theta = fit_theta_likelihood(
            data.groups_data,
            data.obs,
            data.exp,
            data.disp_models,
            data.differential,
            mode=mode,
            theta_x=theta_x,
            theta_step=theta_step,
            theta_tail_sd=theta_tail_sd,
            position_chunk_size=position_chunk_size,
            storage_dtype=storage_dtype,
            log_mu_prior=log_mu_prior,
        )
        return data


class ThetaSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_MEAN_SEGMENTATION,
        log_theta_prior=None,
    ):
        data.theta_segmentation = fit_theta_segmentation(
            data.theta,
            length_prior,
            config,
            log_theta_prior,
        )
        return data
