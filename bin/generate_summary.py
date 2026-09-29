import sys

import numpy as np
import pandas as pd

from footprint_tools.stats.differential_bayesian.api import load_data_results


Z_THRESHOLD = 1.0
PROBABILITY_CUTOFF = 0.99
FOOTPRINT_CUTOFF = 1.0


def posterior_mass(p):
    if hasattr(p, "log_mass"):
        return np.exp(p.log_mass)
    if hasattr(p, "mass"):
        return np.asarray(p.mass)
    raise AttributeError("posterior must expose .log_mass or .mass")


def posterior_grid(p):
    for attr in ("x", "grid", "states", "z_x"):
        if hasattr(p, attr):
            return np.asarray(getattr(p, attr))
    raise AttributeError("posterior must expose a grid")


def prob_gt(p, threshold):
    if hasattr(p, "prob_gt"):
        return np.asarray(p.prob_gt(threshold))
    if hasattr(p, "probability_z_greater_than"):
        return np.asarray(p.probability_z_greater_than(threshold))

    x = posterior_grid(p)
    return posterior_mass(p)[..., x > threshold].sum(axis=-1)


def prob_lt(p, threshold):
    if hasattr(p, "prob_lt"):
        return np.asarray(p.prob_lt(threshold))
    if hasattr(p, "probability_z_less_than"):
        return np.asarray(p.probability_z_less_than(threshold))

    x = posterior_grid(p)
    return posterior_mass(p)[..., x < threshold].sum(axis=-1)


def as_track(x):
    x = np.asarray(x)
    return x[0] if x.ndim == 2 and x.shape[0] == 1 else x


def runs(mask):
    """Return (start, end, value) runs for a boolean 1D array."""
    mask = np.asarray(mask, dtype=bool)

    if not len(mask):
        return []

    change = np.flatnonzero(mask[1:] != mask[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(mask)]

    return list(zip(starts, ends, mask[starts]))


def n_true_runs(mask):
    return sum(value for _, _, value in runs(mask))


def summarize_regions(
    data,
    threshold=Z_THRESHOLD,
    probability_cutoff=PROBABILITY_CUTOFF,
    footprint_cutoff=FOOTPRINT_CUTOFF,
):
    kfp = as_track(data.kfp_zero.mean)
    icc = as_track(data.eta_segmentation.icc.mean)

    posterior = data.common_coefficient_segmentation.posterior
    group_names = tuple(map(str, data.common_coefficient_likelihood.group_names))

    # group x position
    pos = prob_gt(posterior, threshold) >= probability_cutoff
    neg = prob_lt(posterior, -threshold) >= probability_cutoff

    if pos.shape != neg.shape or pos.shape[1] != len(kfp):
        raise ValueError(
            f"shape mismatch: pos={pos.shape}, neg={neg.shape}, kfp={kfp.shape}"
        )

    rows = []

    for region_id, (start, end, is_footprinted) in enumerate(
        runs(kfp >= footprint_cutoff)
    ):
        sl = slice(start, end)

        # Unique bases where >=1 group is significantly positive/negative.
        any_pos = pos[:, sl].any(axis=0)
        any_neg = neg[:, sl].any(axis=0)
        any_diff = any_pos | any_neg

        row = {
            "region_id": region_id,
            "start": start,
            "end": end,
            "width": end - start,
            "status": "footprinted" if is_footprinted else "not_footprinted",
            "kfp_zero_mean": float(kfp[sl].mean()),
            "segmented_icc_mean": float(icc[sl].mean()),

            # Across all groups: unique differential bases.
            "diff_pos_n_bases": int(any_pos.sum()),
            "diff_neg_n_bases": int(any_neg.sum()),
            "diff_n_bases": int(any_diff.sum()),

            # Contiguous differential stretches within this summary region.
            "diff_pos_n_regions": n_true_runs(any_pos),
            "diff_neg_n_regions": n_true_runs(any_neg),
            "diff_n_regions": n_true_runs(any_diff),
        }

        # Per-group information actually needed downstream.
        for i, group in enumerate(group_names):
            row[f"{group}__pos_n_bases"] = int(pos[i, sl].sum())
            row[f"{group}__neg_n_bases"] = int(neg[i, sl].sum())

        rows.append(row)

    return pd.DataFrame(rows)


if __name__ == "__main__":
    dhs_id, npz_path, output_path = sys.argv[1:4]

    data = load_data_results(npz_path)
    summary = summarize_regions(data)

    summary.insert(0, "dhs_id", dhs_id)
    summary.to_csv(output_path, sep="\t", index=False)