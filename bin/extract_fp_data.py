import argparse
from dataclasses import replace

import pandas as pd
from genome_tools import GenomicInterval, df_to_genomic_intervals
from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.basic import GroupsDataLoader
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from footprint_tools.stats.differential_bayesian import DEFAULT_DF_CONFIG
from footprint_tools.stats.differential_bayesian.api import DifferentialFootprintLoader, save_data_results


def load_counts(interval: GenomicInterval, sample_data: pd.DataFrame) -> DataBundle:
    data = DataBundle(interval=interval)
    data = GroupsDataLoader()._load(data, sample_data["extended_annotation"].copy())
    return FootprintsDataLoader()._load(data, footprints_metadata=sample_data, calc_posteriors=False)


def load_dhs(dhs_id: str, dhs_index: str, sample_data: pd.DataFrame) -> DataBundle:
    dhs = pd.read_table(dhs_index).set_index("dhs_id").loc[[dhs_id]].reset_index()
    return load_counts(df_to_genomic_intervals(dhs, extra_columns=["dhs_id"])[0], sample_data)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("dhs_id", "dhs_index", "metadata", "output"):
        p.add_argument(name)
    p.add_argument("--w_min", type=float, default=DEFAULT_DF_CONFIG.w_min, help="call threshold for W")
    a = p.parse_args()
    data = load_dhs(a.dhs_id, a.dhs_index, pd.read_table(a.metadata).set_index("sample_id"))
    save_data_results(DifferentialFootprintLoader()._load(data, replace(DEFAULT_DF_CONFIG, w_min=a.w_min)), a.output)
