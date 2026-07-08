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
from .variance_ratio import VarianceRatioModel


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
