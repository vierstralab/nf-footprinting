from genome_tools import df_to_genomic_intervals

from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from differential.api import DifferentialLoader, GroupMeanSegmentationLoader, VarianceRatioLoader, EtaSegmentationLoader, CoefficientLikelihoodLoader, CoefficientSegmentationLoader
from differential.differential import DifferentialConfig, LengthPrior
from differential.variance_ratio import VarianceRatioConfig
from differential.coefficients import CoefficientConfig


import sys
import pandas as pd
import numpy as np


def extract_group_means(data, length_prior):
    data = DifferentialLoader()._load(
        data,
        config=DifferentialConfig(
            mu_min=-6.0,
            mu_max=6.0,
            n_mu=201,
            n_sig2=81  
        )
    )

    data = GroupMeanSegmentationLoader()._load(
        data,
        length_prior=length_prior,
    )
    return data


def extract_icc(data, length_prior):
    data = VarianceRatioLoader()._load(data, config=VarianceRatioConfig(
        eta_min=-4.0,
        eta_max=4.0,
        eta_step=0.05,
        include_zero=False,
        method="gaussian",
        consistent_mass=None,
        eta_prior_mean=-1.0,
        eta_prior_sd=-2.0,
    ))
    data = EtaSegmentationLoader()._load(
        data,
        length_prior=length_prior,
    )

    return data


def extract_coefs(data, length_prior):
    data = CoefficientLikelihoodLoader()._load(
        data,
        config=CoefficientConfig(
            z_min=-5.0,
            z_max=5.0,
            z_step=0.05,
            method="gaussian",
            zero_mass=None,
            z_prior_sd=1.5,
        )
    )
    data = CoefficientSegmentationLoader()._load(
        data,
        length_prior=length_prior,
    )

    return data

def extract_data_for_dhs_interval(interval, sample_data):
    data = DataBundle(interval=interval)

    data.groups_data = sample_data.loc[:, 'extended_annotation'].copy()
    data.grouping_column = 'extended_annotation'

    data = FootprintsDataLoader()._load(
        data,
        footprints_metadata=sample_data,
        calc_posteriors=False
    )
    print(data.obs.shape)

    length_prior = LengthPrior.from_log_mass(
        np.r_[np.full(5, -np.inf), np.full(40, -1), -1.1],#, np.full(50, -np.inf)],
        infer_tail=True,
    )

    extract_group_means(data, length_prior)
    extract_icc(data, length_prior)
    extract_coefs(data, length_prior)

    return data


if __name__ == "__main__":
    dhs_id = sys.argv[1]  # 'chunk2499_1915_8'
    dhs_index = pd.read_table(sys.argv[2]).set_index('dhs_id')

    sample_data = pd.read_table(sys.argv[3]).set_index('sample_id')

    # for dhs_id in human_anndata.var.query('num_samples >= 4000').index[4:5]:
    dhs_region = df_to_genomic_intervals(
        dhs_index.loc[[dhs_id]].reset_index(),
        extra_columns=['dhs_id']
    )[0]

    data = extract_data_for_dhs_interval(
        dhs_region,
        sample_data
    )
    save_map = {
        'diff_data': data.differential,
        'diff_data_segmentation': data.segmentation,

        'icc_pointwise': data.variance_ratio,
        'icc_segmentation': data.eta_segmentation,

        'coefs_likelihood': data.group_coefficient_likelihood,
        'coefs_segmentation': data.group_coefficient_segmentation,
    }

    for prefix, save_object in save_map.items():
        save_object.to_npz(f'{sys.argv[4]}/{prefix}.{dhs_id}.npz')

