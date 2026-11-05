import argparse

import pandas as pd
from genome_tools import GenomicInterval, df_to_genomic_intervals
from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.basic import GroupsDataLoader
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from footprint_tools.stats.differential_bayesian.api import DifferentialFootprintLoader, save_data_results


def load_counts(interval: GenomicInterval, sample_data: pd.DataFrame) -> DataBundle:
    data = DataBundle(interval=interval)
    data = GroupsDataLoader()._load(data, sample_data["extended_annotation"].copy())
    return FootprintsDataLoader()._load(data, footprints_metadata=sample_data, calc_posteriors=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("dhs_id", "dhs_index", "metadata", "output"):
        p.add_argument(name)
    a = p.parse_args()
    dhs = pd.read_table(a.dhs_index).set_index("dhs_id").loc[[a.dhs_id]].reset_index()
    interval = df_to_genomic_intervals(dhs, extra_columns=["dhs_id"])[0]
    data = load_counts(interval, pd.read_table(a.metadata).set_index("sample_id"))
    save_data_results(DifferentialFootprintLoader()._load(data), a.output)
