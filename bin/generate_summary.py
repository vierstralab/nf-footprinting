import sys
from footprint_tools.stats.differential_bayesian.api import load_data_results

import numpy as np
import pandas as pd

def _posterior_mass(p):
    if hasattr(p, "log_mass"):
        return np.exp(p.log_mass)
    if hasattr(p, "mass"):
        return np.asarray(p.mass)
    raise AttributeError("posterior object must expose `.log_mass` or `.mass`")


def _posterior_grid(p):
    for attr in ("x", "grid", "states", "z_x"):
        if hasattr(p, attr):
            return np.asarray(getattr(p, attr))
    raise AttributeError("posterior object must expose `.x`, `.grid`, `.states`, or `.z_x`")


def _posterior_prob_gt(p, threshold):
    if hasattr(p, "prob_gt"):
        return np.asarray(p.prob_gt(threshold))
    if hasattr(p, "probability_z_greater_than"):
        return np.asarray(p.probability_z_greater_than(threshold))

    x = _posterior_grid(p)
    return _posterior_mass(p)[..., x > float(threshold)].sum(axis=-1)


def _posterior_prob_lt(p, threshold):
    if hasattr(p, "prob_lt"):
        return np.asarray(p.prob_lt(threshold))
    if hasattr(p, "probability_z_less_than"):
        return np.asarray(p.probability_z_less_than(threshold))

    x = _posterior_grid(p)
    return _posterior_mass(p)[..., x < float(threshold)].sum(axis=-1)


def _posterior_prob_abs_gt(p, threshold):
    if hasattr(p, "prob_abs_gt"):
        return np.asarray(p.prob_abs_gt(threshold))
    if hasattr(p, "probability_abs_z_greater_than"):
        return np.asarray(p.probability_abs_z_greater_than(threshold))

    x = _posterior_grid(p)
    return _posterior_mass(p)[..., np.abs(x) > float(threshold)].sum(axis=-1)


def _as_track(x, name):
    x = np.asarray(x)

    if x.ndim == 2 and x.shape[0] == 1:
        x = x[0]

    if x.ndim != 1:
        raise ValueError(
            f"{name} must be 1D or singleton-track x position; got {x.shape}"
        )

    return x


def _region_runs(mask):
    mask = np.asarray(mask, dtype=bool)

    if mask.size == 0:
        return []

    change = np.flatnonzero(mask[1:] != mask[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, mask.size]

    return list(zip(starts, ends, mask[starts]))


def summarize_kfp_zero_regions(
    data,
    *,
    thresholds=(0.85, 1.0, 1.25),
    probability_cutoff=0.99,
    footprint_cutoff=1.0,
    group_names=None,
):
    """
    Region summary from contiguous stretches of data.kfp_zero.mean >= footprint_cutoff.

    Required on data:
        data.kfp_zero
        data.eta_segmentation.icc
        data.common_coefficient_segmentation.posterior

    Regions are 0-based, end-exclusive, relative to the DHS window.

    Main distinction:
        *_n_bases columns count bases.
        n_groups_with_any_* columns count groups with >=1 qualifying base
        inside the region.
    """
    kfp = data.kfp_zero
    icc = data.eta_segmentation.icc
    zpost = data.common_coefficient_segmentation.posterior

    kfp_mean = _as_track(kfp.mean, "kfp_zero.mean")
    icc_mean = _as_track(icc.mean, "eta_segmentation.icc.mean")
    z_mean = np.asarray(zpost.mean)

    if z_mean.ndim != 2:
        raise ValueError(
            f"z posterior mean must be group x position; got {z_mean.shape}"
        )

    n_group, n_pos = z_mean.shape

    if kfp_mean.size != n_pos or icc_mean.size != n_pos:
        raise ValueError(
            f"position mismatch: kfp={kfp_mean.size}, "
            f"icc={icc_mean.size}, z={n_pos}"
        )

    if group_names is None:
        if hasattr(data, "common_coefficient_likelihood"):
            group_names = data.common_coefficient_likelihood.group_names
        elif hasattr(data, "differential"):
            group_names = data.differential.group_names
        else:
            group_names = tuple(f"group_{i}" for i in range(n_group))

    group_names = tuple(map(str, group_names))

    if len(group_names) != n_group:
        raise ValueError(
            f"{len(group_names)=}, but z posterior has {n_group} groups"
        )

    thresholds = tuple(float(x) for x in thresholds)

    p_plus = {
        z_tr: _posterior_prob_gt(zpost, z_tr)
        for z_tr in thresholds
    }

    p_minus = {
        z_tr: _posterior_prob_lt(zpost, -z_tr)
        for z_tr in thresholds
    }

    p_abs = {
        z_tr: _posterior_prob_abs_gt(zpost, z_tr)
        for z_tr in thresholds
    }

    rows = []

    for region_id, (start, end, is_footprinted) in enumerate(
        _region_runs(kfp_mean >= footprint_cutoff)
    ):
        sl = slice(start, end)
        width = int(end - start)

        row = {
            "region_id": int(region_id),
            "start": int(start),
            "end": int(end),
            "width": width,
            "status": "footprinted" if is_footprinted else "not_footprinted",
            "kfp_zero_mean": float(kfp_mean[sl].mean()),
            "segmented_icc_mean": float(icc_mean[sl].mean()),
            "z_mean_min_all_groups": float(z_mean[:, sl].min()),
            "z_mean_max_all_groups": float(z_mean[:, sl].max()),
        }

        for z_tr in thresholds:
            tag = f"z{z_tr:g}"

            pp = p_plus[z_tr][:, sl] >= probability_cutoff
            pm = p_minus[z_tr][:, sl] >= probability_cutoff

            # Directional either:
            # P(z > tr) >= cutoff OR P(z < -tr) >= cutoff.
            pdir = pp | pm

            # True abs posterior:
            # P(|z| > tr) >= cutoff.
            pa = p_abs[z_tr][:, sl] >= probability_cutoff

            # Number of bases in the region where at least one group passes.
            row[f"any_group_p_z_gt_{tag}_n_bases"] = int(pp.any(axis=0).sum())
            row[f"any_group_p_z_lt_minus_{tag}_n_bases"] = int(pm.any(axis=0).sum())
            row[f"any_group_p_directional_abs_z_gt_{tag}_n_bases"] = int(
                pdir.any(axis=0).sum()
            )
            row[f"any_group_p_abs_z_gt_{tag}_n_bases"] = int(pa.any(axis=0).sum())

            # Fractions of bases in the region where at least one group passes.
            row[f"any_group_p_z_gt_{tag}_frac_bases"] = (
                row[f"any_group_p_z_gt_{tag}_n_bases"] / width
            )
            row[f"any_group_p_z_lt_minus_{tag}_frac_bases"] = (
                row[f"any_group_p_z_lt_minus_{tag}_n_bases"] / width
            )
            row[f"any_group_p_directional_abs_z_gt_{tag}_frac_bases"] = (
                row[f"any_group_p_directional_abs_z_gt_{tag}_n_bases"] / width
            )
            row[f"any_group_p_abs_z_gt_{tag}_frac_bases"] = (
                row[f"any_group_p_abs_z_gt_{tag}_n_bases"] / width
            )

            # Number of groups with at least one qualifying base in the region.
            row[f"n_groups_with_any_z_gt_{tag}_base"] = int(pp.any(axis=1).sum())
            row[f"n_groups_with_any_z_lt_minus_{tag}_base"] = int(
                pm.any(axis=1).sum()
            )
            row[f"n_groups_with_any_directional_abs_z_gt_{tag}_base"] = int(
                pdir.any(axis=1).sum()
            )
            row[f"n_groups_with_any_abs_z_gt_{tag}_base"] = int(
                pa.any(axis=1).sum()
            )

            # Fraction of groups with at least one qualifying base.
            row[f"frac_groups_with_any_z_gt_{tag}_base"] = float(
                pp.any(axis=1).mean()
            )
            row[f"frac_groups_with_any_z_lt_minus_{tag}_base"] = float(
                pm.any(axis=1).mean()
            )
            row[f"frac_groups_with_any_directional_abs_z_gt_{tag}_base"] = float(
                pdir.any(axis=1).mean()
            )
            row[f"frac_groups_with_any_abs_z_gt_{tag}_base"] = float(
                pa.any(axis=1).mean()
            )

            for gi, group in enumerate(group_names):
                prefix = f"{group}__{tag}"

                row[f"{prefix}__p_z_gt_n_bases"] = int(pp[gi].sum())
                row[f"{prefix}__p_z_lt_minus_n_bases"] = int(pm[gi].sum())
                row[f"{prefix}__p_directional_abs_z_gt_n_bases"] = int(
                    pdir[gi].sum()
                )
                row[f"{prefix}__p_abs_z_gt_n_bases"] = int(pa[gi].sum())

                row[f"{prefix}__p_z_gt_frac_bases"] = float(pp[gi].mean())
                row[f"{prefix}__p_z_lt_minus_frac_bases"] = float(pm[gi].mean())
                row[f"{prefix}__p_directional_abs_z_gt_frac_bases"] = float(
                    pdir[gi].mean()
                )
                row[f"{prefix}__p_abs_z_gt_frac_bases"] = float(pa[gi].mean())

                row[f"{prefix}__has_any_z_gt_base"] = bool(pp[gi].any())
                row[f"{prefix}__has_any_z_lt_minus_base"] = bool(pm[gi].any())
                row[f"{prefix}__has_any_directional_abs_z_gt_base"] = bool(
                    pdir[gi].any()
                )
                row[f"{prefix}__has_any_abs_z_gt_base"] = bool(pa[gi].any())

        for gi, group in enumerate(group_names):
            vals = z_mean[gi, sl]
            row[f"{group}__z_mean_min"] = float(vals.min())
            row[f"{group}__z_mean_max"] = float(vals.max())

        rows.append(row)

    return pd.DataFrame(rows)


if __name__ == "__main__":
    dhs_id = sys.argv[1]
    npz_path = sys.argv[2]

    data = load_data_results(npz_path)

    summary = summarize_kfp_zero_regions(
        data,
        thresholds=(1.0,),
        probability_cutoff=0.99,
        footprint_cutoff=1.0,
    )

    summary['dhs_id'] = dhs_id
    summary.to_csv(sys.argv[3], index=False, sep='\t')
