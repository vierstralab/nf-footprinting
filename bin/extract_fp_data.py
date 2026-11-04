import argparse

import numpy as np
import pandas as pd
from genome_tools import GenomicInterval, df_to_genomic_intervals
from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.basic import GroupsDataLoader
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from footprint_tools.stats.differential_bayesian import (DEFAULT_DF_MODEL, DifferentialFootprintModelConfig, load_maps,
                                                         profile_stats)
from footprint_tools.stats.differential_bayesian.api import (DifferentialFootprintLoader, DifferentialFootprintSummaryLoader,
                                                             save_data_results)
from footprint_tools.stats.differential_bayesian.hyper import stats


def load_counts(interval: GenomicInterval, sample_data: pd.DataFrame) -> DataBundle:
    data = DataBundle(interval=interval)
    data = GroupsDataLoader()._load(data, sample_data["extended_annotation"].copy())
    return FootprintsDataLoader()._load(data, footprints_metadata=sample_data, calc_posteriors=False)


def profile(data: DataBundle) -> dict[str, np.ndarray]:
    """Binned observed / expected cuts per sample, for fitting the sample maps (samples.py)."""
    g = data.groups_data
    O, E = profile_stats(data.obs.loc[g.index].values, data.exp.loc[g.index].values, g.index.values, g.values)
    return dict(O=O, E=E, sample_ids=np.asarray(g.index, str))


def fit_dhs(data: DataBundle, config: DifferentialFootprintModelConfig = DEFAULT_DF_MODEL,
            maps: pd.DataFrame | None = None) -> tuple[DataBundle, dict[str, np.ndarray]]:
    """Summaries, model fit, and the DHS's statistics for the global constants (hyper.py)."""
    data = DifferentialFootprintSummaryLoader()._load(data, maps=maps)
    data = DifferentialFootprintLoader()._load(data, config=config)
    s = data.summaries
    return data, stats(s.x, s.v, s.group_names, config)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("dhs_id", "dhs_index", "metadata", "output"):
        p.add_argument(name)
    p.add_argument("--profile", action="store_true", help="write only the sample-map profile of the DHS to output")
    p.add_argument("--stats", default="", help="where to write the statistics for hyper")
    p.add_argument("--model_config", default="", help="ModelConfig JSON from hyper (default: package defaults)")
    p.add_argument("--sample_maps", default="", help="TSV of k, f per sample from samples (default: k = 1, f = 0)")
    a = p.parse_args()
    dhs = pd.read_table(a.dhs_index).set_index("dhs_id").loc[[a.dhs_id]].reset_index()
    interval = df_to_genomic_intervals(dhs, extra_columns=["dhs_id"])[0]
    data = load_counts(interval, pd.read_table(a.metadata).set_index("sample_id"))
    if a.profile:
        np.savez_compressed(a.output, **profile(data))
    else:
        config = DifferentialFootprintModelConfig.from_json(a.model_config) if a.model_config else DEFAULT_DF_MODEL
        data, dhs_stats = fit_dhs(data, config, load_maps(a.sample_maps) if a.sample_maps else None)
        save_data_results(data, a.output)
        if a.stats:
            np.savez_compressed(a.stats, **dhs_stats)
