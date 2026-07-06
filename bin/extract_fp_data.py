from genome_tools import df_to_genomic_intervals

from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from differential_test_api import DifferentialFitLoader, WindowedScoreLoader


from collections.abc import Mapping
from pathlib import Path
from typing import Any
import numpy as np
import sys
import pandas as pd


def to_npz(
    data,
    path,
    *,
    compressed: bool=True,
    require_sig2_loglik: bool=True,
    extra: Mapping[str, Any]=None,
) -> Path:
    """Save results needed for plotting and downstream segmentation."""
    path = Path(path)

    if path.suffix != ".npz":
        path = Path(f"{path}.npz")


    fit = data.fit
    likelihood = data.likelihood
    raw_loglik = likelihood.group_loglik_mu_sig2

    if require_sig2_loglik and raw_loglik is None:
        raise ValueError(
            "The position x mu x sig2 likelihood was not retained. "
            "Rerun DifferentialFitLoader with keep_sig2_loglik=True."
        )

    arrays: dict[str, np.ndarray] = {
        # File schema
        "format_version": np.array(1, dtype=np.int16),
        "group_names": np.asarray(
            data.group_names,
            dtype=np.str_,
        ),
        "n_position": np.array(
            fit.statistic.size,
            dtype=np.int64,
        ),

        # Parameter grids
        "mu_x": np.asarray(
            likelihood.mu_x,
            dtype=np.float64,
        ),
        "sig2_x": np.asarray(
            likelihood.sig2_x,
            dtype=np.float64,
        ),

        # Required input for segmentation:
        # group × position × mu
        "group_loglik_mu": np.asarray(
            likelihood.group_loglik_mu,
            dtype=np.float64,
        ),

        # Normalized finite-grid inverse-chi-squared prior masses:
        # group × sig2
        "group_log_sig2_prior_mass": np.asarray(
            likelihood.group_log_sig2_prior_mass,
            dtype=np.float64,
        ),

        # Pointwise differential results
        "lrt_statistic": np.asarray(
            fit.statistic,
            dtype=np.float64,
        ),
        "lrt_df": np.array(
            fit.df,
            dtype=np.int16,
        ),
        "lrt_log_pvalue": np.asarray(
            fit.log_pvalue,
            dtype=np.float64,
        ),
        "common_mu": np.asarray(
            fit.common_mu,
            dtype=np.float64,
        ),
        "group_mu": np.asarray(
            fit.group_mu,
            dtype=np.float64,
        ),
        "group_conditional_sig2": np.asarray(
            fit.group_conditional_sig2,
            dtype=np.float64,
        ),
    }

    # Raw likelihood before applying the sigma² prior:
    # group × position × mu × sig2
    if raw_loglik is not None:
        arrays["group_loglik_mu_sig2"] = np.asarray(
            raw_loglik,
            dtype=np.float64,
        )

    if data.windowed is not None:
        arrays.update({
            "window_radius": np.array(
                data.windowed.radius,
                dtype=np.int32,
            ),
            "window_score": np.asarray(
                data.windowed.score,
                dtype=np.float64,
            ),
            "window_nominal_log_pvalue": np.asarray(
                data.windowed.nominal_log_pvalue,
                dtype=np.float64,
            ),
        })

    if data.diagnostics is not None:
        arrays.update({
            "theta_x": np.asarray(
                data.diagnostics.theta_x,
                dtype=np.float64,
            ),
            "group_invchi2_params": np.asarray(
                data.diagnostics.group_invchi2_params,
                dtype=np.float64,
            ),
            "theta_coverage": np.array(
                data.diagnostics.theta_coverage,
                dtype=np.float64,
            ),
            "mu_boundary_fraction": np.array(
                data.diagnostics.mu_boundary_fraction,
                dtype=np.float64,
            ),
        })

    if extra is not None:
        for key, value in extra.items():
            if not isinstance(key, str) or not key:
                raise ValueError(
                    "Extra metadata keys must be nonempty strings"
                )

            if key in arrays:
                raise KeyError(
                    f"Extra metadata key {key!r} conflicts with "
                    "a reserved result key"
                )

            array = np.asarray(value)

            if array.dtype == object:
                flat = array.ravel()

                if all(
                    isinstance(item, (str, np.str_))
                    for item in flat
                ):
                    array = array.astype(np.str_)
                else:
                    raise TypeError(
                        f"Extra value {key!r} has object dtype. "
                        "Pass a numeric, boolean, or string array instead."
                    )

            arrays[key] = array

    save = np.savez_compressed if compressed else np.savez
    save(path, **arrays)

    return path


def extract_data_for_dhs_interval(interval, sample_data):
    data = DataBundle(interval=interval)

    data.grouping_column = sample_data.loc[:, 'extended_annotation'].copy()

    data = FootprintsDataLoader()._load(
        data,
        footprints_metadata=sample_data,
        calc_posteriors=False
    )
    print(data.obs.shape)
    data = DifferentialFitLoader()._load(data, keep_sig2_loglik=True)
    print('diff loaded')
    data = WindowedScoreLoader()._load(data)

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

    to_npz(data=data, path=sys.argv[4])
    # dhs_region = GenomicInterval("chr15", 74995300, 74995600)



