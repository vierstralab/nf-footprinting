"""Pool per-DHS p-values (diff_pvalues.<dhs_id>.npz): Benjamini-Hochberg over all bases and over DHSs, Storey pi0,
calls, and a page of histograms. python diff_pvalue_report.py FOLDER_OR_FILES... --q 0.05"""
import argparse, glob, os

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from diff_pvalues import EFFECTS, WIN, set_weights

LAM = .5
LABEL = {'any': 'any difference (W = ½)', 'w60': 'W ≥ 0.6', 'w70': 'W ≥ 0.7', 'w80': 'W ≥ 0.8'}


def bh_cutoff(lp, q):
    """log10 of the largest p-value Benjamini-Hochberg rejects at level q."""
    s = np.sort(lp); k = np.flatnonzero(s <= np.log10(q * np.arange(1, len(s) + 1) / len(s)))
    return s[k[-1]] if len(k) else -np.inf


def bh_q(lp):
    """log10 BH q-values."""
    o = np.argsort(lp); n = len(lp); q = lp[o] + np.log10(n / np.arange(1, n + 1))
    out = np.empty(n); out[o] = np.minimum(np.minimum.accumulate(q[::-1])[::-1], 0); return out


def pi0(lp):
    return min(1., (lp > np.log10(LAM)).mean() / (1 - LAM))


def runs(mask):
    e = np.flatnonzero(np.diff(np.r_[0, np.asarray(mask, int), 0])); return zip(e[::2], e[1::2])


def calls(h, d, cut, key='w80'):
    """Runs of >= WIN bases passing the pooled cutoff, split by direction; groups = the set that most often attains
    the statistic over the run; W oriented so that > 0.5 means the set is on the called side."""
    names = [str(x) for x in d['names']]; sets = set_weights(len(names))[0]; more = d[f'w_{key}'] > .5; out = []
    for kind, mask in (('more', more), ('less', ~more)):
        for s, e in runs((d[f'lp_{key}'] <= cut) & mask):
            if e - s < WIN: continue
            S = sets[np.bincount(d[f'set_{key}'][s:e]).argmax()]; w = d[f'w_{key}'][s:e]
            out.append(dict(dhs_id=h, interval=str(d['interval']), kind=kind, start=s, end=e, groups=' + '.join(names[g] for g in S),
                            W=round(float(np.mean(w if kind == 'more' else 1 - w)), 3), min_log10p=round(float(d[f'lp_{key}'][s:e].min()), 2)))
    return out


def figure(D, base, dhs5, cut, q, path):
    fig, axs = plt.subplots(2, 3, figsize=(15, 9.5), gridspec_kw=dict(hspace=.5, wspace=.28, left=.06, right=.97, top=.88, bottom=.08))
    fig.suptitle(f'Differential footprint p-values: {len(D):,} DHSs, {len(base["any"]):,} bases; BH at {q:g}', x=.06, ha='left', fontsize=12, weight='bold')
    bins = np.linspace(0, 1, 21)
    for ax, k in zip(axs[0, :2], ('any', 'w80')):
        lp = base[k]; ax.hist(10 ** lp, bins=bins, color='0.6', density=True); e = pi0(lp)
        ax.axhline(e, c='k', ls='--', lw=1); ax.axvline(LAM, c='0.4', ls=':', lw=.8)
        ax.set_title(f'per base · null: {LABEL[k]} not met\nStorey π0 = {e:.2f}; BH rejects {(lp <= cut[k]).mean():.1%} of bases', fontsize=9); ax.set_xlabel('p')
    ax = axs[0, 2]
    for k, c in (('any', '0.5'), ('w80', '#1f5bff')):
        lp = dhs5[k]; ax.hist(10 ** lp, bins=bins, histtype='step', lw=1.6, color=c, density=True,
                              label=f'{LABEL[k]}: π0 {pi0(lp):.2f}, BH rejects {(bh_q(lp) <= np.log10(q)).mean():.1%}')
    ax.legend(fontsize=7.5, frameon=False, loc='upper center'); ax.set_xlabel('p'); ax.set_title('per DHS (best 5-bp window)', fontsize=9)
    ax = axs[1, 0]; x = list(EFFECTS)
    ax.plot(x, [(bh_q(dhs5[k]) <= np.log10(q)).mean() for k in x], 'o-', c='k', label='DHSs: BH discoveries')
    ax.plot(x, [1 - pi0(dhs5[k]) for k in x], 'o:', c='k', label='DHSs: 1 − π0')
    ax.plot(x, [(base[k] <= cut[k]).mean() for k in x], 's-', c='0.55', label='bases: BH discoveries')
    ax.plot(x, [1 - pi0(base[k]) for k in x], 's:', c='0.55', label='bases: 1 − π0')
    ax.set_xticks(range(len(x)), [LABEL[k] for k in x], fontsize=7.5); ax.set_ylim(0, 1.02); ax.legend(fontsize=7.5, frameon=False)
    ax.set_title('Fraction differential vs the effect required\n(above "any", 1 − π0 is a lower bound: nulls below the threshold give p near 1)', fontsize=9)
    ax = axs[1, 1]; xi = np.concatenate([d['xi'] for d in D.values()]); ax.hist(np.clip(xi, -1, 1), bins=50, color='0.6'); ax.axvline(0, c='k', lw=.8)
    ax.set_title(f'GPD shape ξ of the permutation maxima, per base\n{(xi > 0).mean():.0%} > 0 (heavier than exponential); ξ < 0 is used as 0', fontsize=9); ax.set_xlabel('ξ')
    ax = axs[1, 2]; c95 = np.concatenate([d['c95'] for d in D.values()]); ax.hist(np.clip(c95, 0, 8), bins=50, color='0.6')
    ax.set_title(f'95% point of the permutation maximum, per base: median {np.median(c95):.2f}', fontsize=9); ax.set_xlabel('critical value')
    fig.savefig(path)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('inputs', nargs='+'); p.add_argument('--q', type=float, default=.05)
    p.add_argument('--prefix', default='diff_pvalues')
    a = p.parse_args()
    files = sorted(f for x in a.inputs for f in (glob.glob(f'{x}/diff_pvalues.*.npz') if os.path.isdir(x) else [x]))
    D = {os.path.basename(f)[len('diff_pvalues.'):-4]: dict(np.load(f)) for f in files}
    base = {k: np.concatenate([d[f'lp_{k}'] for d in D.values()]) for k in EFFECTS}
    dhs5 = {k: np.array([float(d[f'dhs5_{k}']) for d in D.values()]) for k in EFFECTS}
    cut = {k: bh_cutoff(base[k], a.q) for k in EFFECTS}
    table = pd.DataFrame({'dhs_id': list(D), 'interval': [str(d['interval']) for d in D.values()]})
    for k in EFFECTS:
        table[f'log10p_{k}'], table[f'log10q_{k}'] = dhs5[k].round(3), bh_q(dhs5[k]).round(3)
        table[f'n_bases_bh_{k}'] = [int((d[f'lp_{k}'] <= cut[k]).sum()) for d in D.values()]
    table.to_csv(f'{a.prefix}.dhs.tsv', sep='\t', index=False)
    pd.DataFrame([c for h, d in D.items() for c in calls(h, d, cut['w80'])]).to_csv(f'{a.prefix}.calls_w80.tsv', sep='\t', index=False)
    figure(D, base, dhs5, cut, a.q, f'{a.prefix}.summary.pdf')
    with open(f'{a.prefix}.cutoffs.tsv', 'w') as fh:
        fh.write('effect\tbase_log10p_cutoff\tbase_frac_bh\tbase_pi0\tdhs_frac_bh\tdhs_pi0\n')
        for k in EFFECTS:
            fh.write(f'{k}\t{cut[k]:.3f}\t{(base[k] <= cut[k]).mean():.4f}\t{pi0(base[k]):.3f}\t{(bh_q(dhs5[k]) <= np.log10(a.q)).mean():.4f}\t{pi0(dhs5[k]):.3f}\n')
