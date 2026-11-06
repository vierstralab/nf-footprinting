"""Diagnostic report for one DHS from its saved differential footprint fit: PDF, candidate table, one-row DHS summary."""
import argparse

import numpy as np
import pandas as pd

from footprint_tools.stats.differential_bayesian.api import load_data_results
from diff_debug.diag import DHS, candidate_table, dhs_summary
from diff_debug.plots import report
from extract_fp_data import load_dhs

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("dhs_id", "dhs_index", "metadata", "results", "pdf", "candidates", "summary"):
        p.add_argument(name)
    p.add_argument("--rule", choices=("base", "window"), default="base", help="call rule the report treats as the calls")
    p.add_argument("--spot_column", default="SPOT3_score")
    a = p.parse_args()
    meta = pd.read_table(a.metadata).set_index("sample_id")
    data = load_dhs(a.dhs_id, a.dhs_index, meta)
    F = load_data_results(a.results, data).footprints
    s = list(F.sample_names)
    spot = meta[a.spot_column].reindex(s).to_numpy(float) if a.spot_column in meta else np.full(len(s), np.nan)
    d = DHS(F, data.obs.loc[s].to_numpy(float), data.exp.loc[s].to_numpy(float), list(data.disp_models.loc[s]),
            spot, a.dhs_id, data.interval.to_ucsc(), a.rule)
    C = report(d, a.pdf)
    candidate_table(d, C).to_csv(a.candidates, sep="\t", index=False)
    dhs_summary(d, C).to_csv(a.summary, sep="\t", index=False)
