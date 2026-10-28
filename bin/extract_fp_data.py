import sys

import pandas as pd
from genome_tools import df_to_genomic_intervals
from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.basic import GroupsDataLoader
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from footprint_tools.stats.differential_bayesian import DEFAULT_NESTED, NestedConfig
from footprint_tools.stats.differential_bayesian.api import NestedLoader, SummaryLoader, save_data_results


def fit_dhs(interval, sample_data, config=DEFAULT_NESTED):
    data = DataBundle(interval=interval)
    data = GroupsDataLoader()._load(data, sample_data["extended_annotation"].copy())
    data = FootprintsDataLoader()._load(data, footprints_metadata=sample_data, calc_posteriors=False)
    data = SummaryLoader()._load(data)
    return NestedLoader()._load(data, config=config)


if __name__ == "__main__":
    # dhs_id dhs_index sample_metadata output.npz [nested_config.json]
    dhs_id, dhs_index, metadata, output = sys.argv[1:5]
    config = NestedConfig.from_json(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] else DEFAULT_NESTED
    dhs = pd.read_table(dhs_index).set_index("dhs_id").loc[[dhs_id]].reset_index()
    interval = df_to_genomic_intervals(dhs, extra_columns=["dhs_id"])[0]
    save_data_results(fit_dhs(interval, pd.read_table(metadata).set_index("sample_id"), config), output)
