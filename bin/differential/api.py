from genome_tools.plotting.modular_plot.api import PlotDataLoader

from .coefficients import (
    CommonCoefficientModel,
    ZeroCoefficientModel,
    fit_coefficient_segmentation,
    infer_kfp_dev,
    infer_kfp_zero,
)
from .config import (
    DEFAULT_COEFFICIENT,
    DEFAULT_COEFFICIENT_SEGMENTATION,
    DEFAULT_DIFFERENTIAL,
    DEFAULT_ETA_SEGMENTATION,
    DEFAULT_MEAN_SEGMENTATION,
    DEFAULT_THETA,
    DEFAULT_THETA_SEGMENTATION,
    DEFAULT_VARIANCE_RATIO,
)
from .differential import DifferentialModel, fit_group_means_segmentation
from .eta import fit_eta_segmentation
from .theta import ThetaModel, fit_theta_segmentation
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


class CommonCoefficientLikelihoodLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_COEFFICIENT,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.common_coefficient_likelihood = CommonCoefficientModel(config).fit(
            data.differential,
            data.variance_ratio,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class CommonCoefficientSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_COEFFICIENT_SEGMENTATION,
    ):
        data.common_coefficient_segmentation = fit_coefficient_segmentation(
            data.common_coefficient_likelihood,
            length_prior,
            config,
        )
        return data


class ZeroCoefficientLikelihoodLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_COEFFICIENT,
        log_z_prior=None,
    ):
        data.zero_coefficient_likelihood = ZeroCoefficientModel(config).fit(
            data.differential,
            log_z_prior,
        )
        return data


class ZeroCoefficientSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_COEFFICIENT_SEGMENTATION,
    ):
        data.zero_coefficient_segmentation = fit_coefficient_segmentation(
            data.zero_coefficient_likelihood,
            length_prior,
            config,
        )
        return data


class ZeroFootprintCountLoader(PlotDataLoader):
    def _load(self, data, threshold=1.0, source="segmented"):
        if source == "segmented":
            posterior = data.zero_coefficient_segmentation.posterior
        elif source == "pointwise":
            posterior = data.zero_coefficient_likelihood.posterior()
        else:
            raise ValueError("source must be 'segmented' or 'pointwise'")
        data.kfp_zero = infer_kfp_zero(posterior, threshold)
        return data


class DeviantFootprintCountLoader(PlotDataLoader):
    def _load(self, data, threshold=1.0, source="segmented"):
        if source == "segmented":
            posterior = data.common_coefficient_segmentation.posterior
        elif source == "pointwise":
            posterior = data.common_coefficient_likelihood.posterior()
        else:
            raise ValueError("source must be 'segmented' or 'pointwise'")
        data.kfp_dev = infer_kfp_dev(posterior, threshold)
        return data


class ThetaLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_THETA,
        log_mu_prior=None,
    ):
        data.theta = ThetaModel(config).fit(
            data.groups_data,
            data.obs,
            data.exp,
            data.disp_models,
            data.differential,
            log_mu_prior,
        )
        return data


class ThetaSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_THETA_SEGMENTATION,
        log_theta_prior=None,
    ):
        data.theta_segmentation = fit_theta_segmentation(
            data.theta,
            length_prior,
            config,
            log_theta_prior,
        )
        return data
