from genome_tools import df_to_genomic_intervals

from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from differential.api import (
    DifferentialLoader, GroupMeanSegmentationLoader, VarianceRatioLoader, EtaSegmentationLoader,
    ZeroCoefficientLikelihoodLoader, ZeroCoefficientSegmentationLoader,
    CommonCoefficientLikelihoodLoader, CommonCoefficientSegmentationLoader,
    ZeroFootprintCountLoader, ThetaLoader, Mu0SegmentationLoader
)
from differential.segmentation import LengthPrior
from differential.config import VarianceRatioConfig, DifferentialConfig, CoefficientConfig, ThetaConfig

from differential_stats import summarize_kfp_zero_regions


import sys
import pandas as pd
import numpy as np
import polars as pl


def extract_group_means(data, length_prior, config):
    data = DifferentialLoader()._load(data, config=config)
    data = GroupMeanSegmentationLoader()._load(data, length_prior=length_prior)

    return data


def extract_per_sample_depletions(data, config):
    data = ThetaLoader()._load(data, config=config)

    return data


def extract_icc(data, length_prior, config):
    data = VarianceRatioLoader()._load(data, config=config)
    data = Mu0SegmentationLoader()._load(data, length_prior=length_prior)
    data = EtaSegmentationLoader()._load(data, length_prior=length_prior)

    return data


def extract_zero_coefs(data, length_prior, config):
    data = ZeroCoefficientLikelihoodLoader()._load(data, config=config)
    data = ZeroCoefficientSegmentationLoader()._load(data, length_prior=length_prior)

    return data

def extract_common_coefs(data, length_prior, config):
    data = CommonCoefficientLikelihoodLoader()._load(data, config=config)
    data = CommonCoefficientSegmentationLoader()._load(data, length_prior=length_prior)

    return data

def footprint_log_mass(path, max_width=100, min_width=4):
    df = pl.read_csv(path, separator="\t")
    hs = ["hotspot_chr", "hotspot_start", "hotspot_end"]

    footprint_widths = (
        df.unique("dhs_id")
        .select(pl.col("dhs_width").alias("width"))
    )

    spacings = (
        df.unique(hs + ["dhs_id"])
        .sort(hs + ["start", "end"])
        .with_columns(
            prev_id=pl.col("dhs_id").shift().over(hs),
            width=(
                pl.col("start") - pl.col("end").shift().over(hs)
            ).clip(lower_bound=0),
        )
        .drop_nulls("prev_id")
        .unique(["prev_id", "dhs_id"])
        .select("width")
    )

    counts_df = (
        pl.concat([footprint_widths, spacings])
        .filter(pl.col("width").is_between(0, max_width))
        .group_by("width")
        .len()
    )

    counts = np.zeros(max_width + 1)
    counts[counts_df["width"].to_numpy().astype(int)] = counts_df["len"].to_numpy()

    log_mass = np.full(max_width + 1, -np.inf)
    positive = counts > 0
    log_mass[positive] = np.log(counts[positive] / counts.sum())
    log_mass[:min_width] = -np.inf

    return log_mass

def extract_data_for_dhs_interval(interval, sample_data, fp_index_w_hotspots_path):
    data = DataBundle(interval=interval)

    data.groups_data = sample_data.loc[:, 'extended_annotation'].copy()
    data.grouping_column = 'extended_annotation'

    data = FootprintsDataLoader()._load(
        data,
        footprints_metadata=sample_data,
        calc_posteriors=False
    )

    length_prior = LengthPrior.from_log_mass(
        footprint_log_mass(
            fp_index_w_hotspots_path,
            max_width=500,
            min_width=4,
        ),
        infer_tail=True,
        n_tail=2,
    )

    data = extract_group_means(
        data, length_prior, config=DifferentialConfig(
            mu_min=-6.0,
            mu_max=6.0,
            n_mu=201,
            n_sig2=81  
        )
    )
    data = extract_per_sample_depletions(
        data, config=ThetaConfig(
            mode="sample_only",
            position_chunk_size=64,
        )
    )
    data = extract_icc(
        data, length_prior, config=VarianceRatioConfig(
            eta_min=-4.0,
            eta_max=4.0,
            eta_step=0.05,
            include_zero=False,
            consistent_mass=None,
            eta_prior_mean=-2.0,
            eta_prior_sd=4.0,
            method="gaussian",
        )
    )
    
    coef_config = CoefficientConfig(
        z_min=-5.0,
        z_max=5.0,
        z_step=0.05,
    )
    data = extract_zero_coefs(
        data, length_prior, config=coef_config,
    )
    data = extract_common_coefs(
        data, length_prior, config=coef_config,
    )

    data = ZeroFootprintCountLoader()._load(
        data,
        threshold=1.0,
        source="segmented",
    )

    return data


if __name__ == "__main__":
    dhs_id = sys.argv[1]  # 'chunk2499_1915_8'
    dhs_index = pd.read_table(sys.argv[2]).set_index('dhs_id')

    sample_data = pd.read_table(sys.argv[3]).set_index('sample_id')

    fp_index_w_hotspots_path = sys.argv[4]

    # for dhs_id in human_anndata.var.query('num_samples >= 4000').index[4:5]:
    dhs_region = df_to_genomic_intervals(
        dhs_index.loc[[dhs_id]].reset_index(),
        extra_columns=['dhs_id']
    )[0]

    data = extract_data_for_dhs_interval(
        dhs_region,
        sample_data,
        fp_index_w_hotspots_path
    )

    save_map = {
        'diff_data': data.differential,
        'diff_data_segmentation': data.segmentation,

        'icc_pointwise': data.variance_ratio,
        'icc_segmentation': data.eta_segmentation,

        'zero_coefs_likelihood': data.zero_coefficient_likelihood,
        'zero_coefs_segmentation': data.zero_coefficient_segmentation,

        'common_coefs_likelihood': data.common_coefficient_likelihood,
        'common_coefs_segmentation': data.common_coefficient_segmentation,

        'per_sample_depletions': data.theta,
        'mu0_segmentation': data.mu0_segmentation,
    }

    summary = summarize_kfp_zero_regions(
        data,
        thresholds=(0.85, 1.0, 1.25),
        probability_cutoff=0.99,
        footprint_cutoff=1.0,
    )

    summary['dhs_id'] = dhs_id
    summary.to_csv(f'{sys.argv[4]}/summary.{dhs_id}.tsv', index=False, sep='\t')

    for prefix, save_object in save_map.items():
        save_object.to_npz(f'{sys.argv[4]}/{prefix}.{dhs_id}.npz')

