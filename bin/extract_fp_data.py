import argparse
import json

import numpy as np
import pandas as pd
import pysam
from genome_tools import df_to_genomic_intervals
from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.basic import GroupsDataLoader
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from footprint_tools.stats.differential_bayesian import DEFAULT_MODEL, ModelConfig, seqcorr
from footprint_tools.stats.differential_bayesian.api import FootprintLoader, SummaryLoader, save_data_results
from footprint_tools.stats.differential_bayesian.hyper import evidence_tables


def fit_dhs(interval, sample_data, fasta, config=DEFAULT_MODEL, coef=None):
    """Fit one DHS; also returns the statistics for the global fits (sequence correction, constants)."""
    data = DataBundle(interval=interval)
    data = GroupsDataLoader()._load(data, sample_data["extended_annotation"].copy())
    data = FootprintsDataLoader()._load(data, footprints_metadata=sample_data, calc_posteriors=False)
    seq = fasta.fetch(interval.chrom, interval.start - seqcorr.W, interval.end + seqcorr.W + 1)
    stats = {"chrom": interval.chrom, "seqcorr.obs": data.obs.values.sum(0), "seqcorr.exp": data.exp.values.sum(0),
             "seqcorr.seq": seq}
    if coef is not None:
        data.exp = data.exp * 2 ** seqcorr.correction(seq, coef)
    data = SummaryLoader()._load(data)
    data = FootprintLoader()._load(data, config=config)
    stats["evidence.consensus"], stats["evidence.pattern"] = evidence_tables(data.summaries.x, data.summaries.v, config)
    return data, stats


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("dhs_id", "dhs_index", "metadata", "fasta", "output", "stats"):
        p.add_argument(name)
    p.add_argument("--model_config", default="", help="ModelConfig JSON from hyper (default: package defaults)")
    p.add_argument("--seqcorr", default="", help="sequence correction JSON from seqcorr (default: none)")
    a = p.parse_args()
    config = ModelConfig.from_json(a.model_config) if a.model_config else DEFAULT_MODEL
    coef = json.load(open(a.seqcorr)) if a.seqcorr else None
    dhs = pd.read_table(a.dhs_index).set_index("dhs_id").loc[[a.dhs_id]].reset_index()
    interval = df_to_genomic_intervals(dhs, extra_columns=["dhs_id"])[0]
    data, stats = fit_dhs(interval, pd.read_table(a.metadata).set_index("sample_id"), pysam.FastaFile(a.fasta), config, coef)
    save_data_results(data, a.output)
    np.savez_compressed(a.stats, **stats)
