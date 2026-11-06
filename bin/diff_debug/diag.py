"""Diagnostics for one DHS fit: candidates (called and near-called), why each is or isn't called, and per-sample QC."""
from dataclasses import dataclass
from functools import cached_property

import numpy as np
import pandas as pd
from scipy.special import ndtr
from scipy.stats import nbinom, theilslopes

from footprint_tools.stats.differential_bayesian.calls import runs
from footprint_tools.stats.differential_bayesian.results import DifferentialFootprintFit
from footprint_tools.stats.differential_bayesian.segmentation import _pwl

from .window import SIGN, base_pass, called_window, cover, window_closeness, window_tracks

RULES = {'base': 'per-base rule', 'window': 'window rule'}


@dataclass
class DHS:
    fit: DifferentialFootprintFit
    obs: np.ndarray            # samples x positions, rows in fit.sample_names order
    exp: np.ndarray
    models: list               # dispersion models, same order
    spot: np.ndarray
    name: str
    interval: str
    rule: str = 'base'         # which call rule the report treats as the calls: 'base' or 'window'

    @cached_property
    def idx(self):
        lab = np.array(self.fit.sample_groups)
        return [np.flatnonzero(lab == g) for g in self.fit.group_names]

    @cached_property
    def others(self):
        """G x n: median of the other groups' means."""
        M = self.fit.mean
        return np.stack([np.median(np.delete(M, g, 0), 0) for g in range(len(M))])

    @cached_property
    def wt(self):
        return window_tracks(self.fit, self.idx)

    @cached_property
    def qc(self):
        return sample_qc(self)

    def wclose(self, kind):
        """Per 5-bp window (index = first base): min(lower bound / a_min, signed Delta / fold) on window-averaged depths."""
        return window_closeness(self.fit, self.wt, kind)

    def called(self, kind, rule=None):
        return called_window(self.fit, self.wt, kind) if (rule or self.rule) == 'window' else self.fit.called(kind)

    @property
    def alt(self):
        return 'base' if self.rule == 'window' else 'window'


def scaling(d: DHS, g: int, exclude: np.ndarray | None = None):
    """Group mean vs the other groups' median over footprinted bases outside `exclude`: robust slope, intercept."""
    keep = d.fit.footprinted & (True if exclude is None else ~exclude)
    return np.array(theilslopes(d.fit.mean[g, keep], d.others[g, keep])[:2]) if keep.sum() > 10 else np.array([1., 0.])


def candidate_mask(d: DHS, near: float = 0.85) -> np.ndarray:
    """G x n: bases in a window within `near` of passing, either kind (excluded when estimating a group's scaling)."""
    return cover(np.maximum(d.wclose('more'), d.wclose('less')) >= near, d.fit.p_any.size, d.wt['w'])


def pairwise(d: DHS, g: int, kind: str, cols: np.ndarray) -> pd.Series:
    """pi_gh on depths averaged over `cols` (as in the window test); high = g is on the called side of h."""
    m, v, idx = d.fit.theta_mean[:, cols].mean(1), d.fit.theta_var[:, cols].mean(1), d.idx
    out = {}
    for h, ih in enumerate(idx):
        if h == g: continue
        P = ndtr((m[ih][None] - m[idx[g]][:, None]) / np.sqrt(v[idx[g]][:, None] + v[ih][None])).mean()
        out[d.fit.group_names[h]] = P if kind == 'more' else 1 - P
    return pd.Series(out).sort_values()


def sample_side(d: DHS, g: int, kind: str, cols: np.ndarray) -> np.ndarray:
    """For each sample of g: P(it is on the called side of a random sample of another group), averaged over groups."""
    m, v, idx = d.fit.theta_mean[:, cols].mean(1), d.fit.theta_var[:, cols].mean(1), d.idx
    return np.mean([ndtr(-SIGN[kind] * (m[ih][None] - m[idx[g]][:, None]) / np.sqrt(v[idx[g]][:, None] + v[ih][None])).mean(1)
                    for h, ih in enumerate(idx) if h != g], 0)


def candidates(d: DHS, near: float = 0.85) -> pd.DataFrame:
    """Each called run, and each run (>= 5 bp, outside calls) covered by windows within `near` of passing.
    Every candidate is reported on its best 5-bp window."""
    F, c = d.fit, d.fit.cfg
    n, w = F.p_any.size, d.wt['w']
    rows, excl = [], candidate_mask(d, near)
    for kind in ('more', 'less'):
        wc, called, alt, base = d.wclose(kind), d.called(kind), d.called(kind, d.alt), base_pass(F, kind)
        A, lo, sd = getattr(F, f'a_{kind}'), getattr(F, f'lo_{kind}'), SIGN[kind] * F.delta
        for g, name in enumerate(F.group_names):
            for s, e in runs(called[g]) + runs(cover(wc[g:g + 1] >= near, n, w)[0] & ~called[g], w):
                s0 = min(s, n - w); ts = np.arange(s0, max(s0, e - w) + 1); t = ts[np.argmax(wc[g, ts])]; win = np.arange(t, t + w)
                r = dict(group=name, g=g, kind=kind, start=s, end=e, called=called[g, s:e].any(), n_called=int(called[g, s:e].sum()),
                         alt_bp=int(alt[g, s:e].sum()), n_pass=int(base[g, s:e].sum()), n_pass_win=int(base[g, win].sum()),
                         longest_pass=max((j - i for i, j in runs(base[g, s:e])), default=0),
                         win=win, a=d.wt[f'a_{kind}'][g, t], lo=d.wt[f'lo_{kind}'][g, t], delta=d.wt['delta'][g, t], close=wc[g, t])
                r['outcome'], r['detail'] = outcome(r, A[g, win], lo[g, win], sd[g, win], c, d.rule)
                pi, q = pairwise(d, g, kind, win), sample_side(d, g, kind, win)
                r['ties'] = ', '.join(f'{k} {v:.2f}' for k, v in pi[pi < c.a_min].head(3).items()) or '-'
                r['n_side'], r['n_g'] = int((q > .5).sum()), len(q)
                slope, icpt = scaling(d, g, excl[g])
                r['slope'], r['delta_scale'] = float(slope), float(np.mean((slope - 1) * d.others[g, win] + icpt))
                ig, io = d.idx[g], np.setdiff1d(np.arange(len(F.sample_names)), d.idx[g])
                r['exp_g'], r['exp_o'] = float(np.median(d.exp[ig][:, win])), float(np.median(d.exp[io][:, win]))
                r['sd_g'] = float(np.median(np.sqrt(F.theta_var[ig][:, win])))
                r['spread_g'] = float(np.median(F.theta_mean[ig][:, win].std(0)))
                r['warn'] = flags(r)
                rows.append(r)
    df = pd.DataFrame(rows)
    if df.empty: return df
    df = df.assign(key=np.where(df.called, 10 + df.lo, df.close)).sort_values('key', ascending=False).drop(columns='key').reset_index(drop=True)
    df['tag'] = [f'C{i + 1}' for i in range(len(df))]
    return df


def outcome(r: dict, a, lb, sd, c, rule) -> tuple[str, list[str]]:
    """Why the candidate is not called (window-level failures, or the per-base run length), and each window base on its own."""
    detail = []
    for j, b in enumerate(r['win']):
        why = ([f'A {a[j]:.3f} < {c.a_min}'] if a[j] < c.a_min else [f'lower bound {lb[j]:.3f} (A {a[j]:.2f})'] if lb[j] < c.a_min else []) + \
              ([f'|Δ| {abs(sd[j]):.3f} ≤ {c.fold:.3f}' if sd[j] > 0 else 'Δ in the wrong direction'] if sd[j] <= c.fold else [])
        detail.append(f"{b}: " + ('; '.join(why) if why else 'passes'))
    if r['called']: return 'called', detail
    parts = (['overlap' if r['a'] < c.a_min else 'uncertain'] if r['lo'] < c.a_min else []) + \
            (['small effect'] if SIGN[r['kind']] * r['delta'] <= c.fold else [])
    if not parts and rule == 'base':
        parts = [f"too short: {r['longest_pass']} bp pass in a row"]
    return ' + '.join(parts), detail


def flags(r: dict) -> list[str]:
    """Warnings worth a look, whether or not the candidate is called."""
    out = []
    if np.sign(r['delta_scale']) == np.sign(r['delta']) and abs(r['delta_scale']) >= .5 * abs(r['delta']):
        out.append('depth scaling explains ≥ ½ of Δ')
    if r['n_side'] < .75 * r['n_g']: out.append(f"only {r['n_side']}/{r['n_g']} samples on the called side")
    if r['exp_g'] < 2: out.append('low coverage')
    if r['sd_g'] > abs(r['delta']): out.append('count noise ≥ effect')
    return out


def nb_params(models, e):
    mp, rp = (np.stack([getattr(m, k) for m in models]) for k in ('mu_params', 'r_params'))
    mu = np.maximum(np.where(e < mp[:, 0, None], mp[:, 6, None] * e, _pwl(e, mp, 3)), 1e-4)
    r = 1 / np.maximum(_pwl(np.maximum(e, rp[:, 0, None]), rp, 5), 0.01)
    return mu, r


def sample_qc(d: DHS, seed: int = 0) -> tuple[pd.DataFrame, dict]:
    """Per-sample checks of the stage-1 assumptions, and pooled PIT values (at the posterior mean depth)."""
    F = d.fit; m = F.theta_mean.astype(float); fp = F.footprinted
    mu, r = nb_params(d.models, d.exp * np.exp2(m))
    y, ok = d.obs, d.exp > 0
    resid2 = (y - mu) ** 2 / (mu + mu ** 2 / r)
    p = r / (r + mu)
    pit = nbinom.cdf(y - 1, r, p) + np.random.default_rng(seed).uniform(size=y.shape) * nbinom.pmf(y, r, p)
    loo = np.empty(len(m))
    for i in d.idx:
        tot = m[i].sum(0)
        for s in i: loo[s] = np.corrcoef(m[s], (tot - m[s]) / (len(i) - 1))[0, 1]
    qc = pd.DataFrame(dict(sample=F.sample_names, group=F.sample_groups, spot=d.spot, cuts=d.obs.sum(1),
                           disp=[resid2[s][ok[s]].mean() for s in range(len(m))],
                           offset=np.median(m[:, ~fp], 1) if (~fp).any() else 0., loo=loo,
                           fp_depth=m[:, fp].mean(1) if fp.any() else np.nan))
    return qc, dict(pit_fp=pit[:, fp][ok[:, fp]], pit_open=pit[:, ~fp][ok[:, ~fp]])


def dhs_flags(d: DHS) -> list[str]:
    """DHS-level warnings from the QC checks."""
    qc, out = d.qc[0], []
    off, disp = qc.offset.median(), qc.disp.median()
    if abs(off) > .5: out.append(f'open DNA sits at {off:+.2f}, not 0: footprints use an absolute threshold')
    if disp > 1.5: out.append(f'counts noisier than the dispersion model (index {disp:.2f})')
    med = qc.groupby('group').loo.transform('median')
    odd = qc[(qc.loo < med - .25) | ((qc.offset - off).abs() > 1)]
    if len(odd): out.append('outlier samples: ' + ', '.join(f'{r.sample} ({r.group})' for r in odd.head(6).itertuples()))
    excl = candidate_mask(d)
    sl = {g: scaling(d, k, excl[k])[0] for k, g in enumerate(d.fit.group_names)}
    off_sl = [f'{g} {v:.2f}' for g, v in sl.items() if abs(v - 1) > .25]
    if off_sl: out.append('group-wide depth scaling: ' + ', '.join(off_sl))
    return out


# ---------- tables ----------

CANDIDATE_COLUMNS = ['dhs_id', 'rule', 'tag', 'group', 'kind', 'start', 'end', 'called', 'called_bp', 'alt_rule_bp', 'outcome',
                     'win_start', 'win_end', 'closeness', 'a', 'lower_bound', 'delta', 'pass_in_window', 'pass_in_run',
                     'ties', 'side', 'n_group', 'scaling_slope', 'scaling_delta', 'exp_group', 'exp_others',
                     'posterior_sd', 'between_sample_sd', 'warnings']


def candidate_table(d: DHS, C: pd.DataFrame) -> pd.DataFrame:
    if C.empty: return pd.DataFrame(columns=CANDIDATE_COLUMNS)
    return pd.DataFrame(dict(
        dhs_id=d.name, rule=d.rule, tag=C.tag, group=C.group, kind=C.kind, start=C.start, end=C.end, called=C.called,
        called_bp=C.n_called, alt_rule_bp=C.alt_bp, outcome=C.outcome, win_start=C.win.map(min), win_end=C.win.map(max) + 1,
        closeness=C.close, a=C.a, lower_bound=C.lo, delta=C.delta, pass_in_window=C.n_pass_win, pass_in_run=C.n_pass,
        ties=C.ties, side=C.n_side, n_group=C.n_g, scaling_slope=C.slope, scaling_delta=C.delta_scale, exp_group=C.exp_g,
        exp_others=C.exp_o, posterior_sd=C.sd_g, between_sample_sd=C.spread_g, warnings=C.warn.map('; '.join)))[CANDIDATE_COLUMNS].round(4)


def dhs_summary(d: DHS, C: pd.DataFrame) -> pd.DataFrame:
    qc = d.qc[0]
    n_calls = lambda rule: sum(len(runs(d.called(k, rule)[g])) for k in ('more', 'less') for g in range(len(d.fit.group_names)))
    return pd.DataFrame([dict(
        dhs_id=d.name, interval=d.interval, rule=d.rule, width=d.fit.p_any.size, n_footprints=len(runs(d.fit.footprinted)),
        n_calls=n_calls(d.rule), n_calls_alt_rule=n_calls(d.alt), n_near=int((~C.called).sum()) if len(C) else 0,
        median_offset=qc.offset.median(), median_dispersion_index=qc.disp.median(), warnings='; '.join(dhs_flags(d)))]).round(4)
