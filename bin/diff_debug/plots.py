"""Report pages for one DHS: overview, candidate table, one row per called / near-called candidate, QC of the assumptions."""
import textwrap
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle

from footprint_tools.stats.differential_bayesian.calls import runs

from .diag import RULES, SIGN, candidate_mask, candidates, dhs_flags, pairwise, sample_side, scaling
from .window import base_pass

CALL = {'more': '#1f5bff', 'less': '#ff2d2d'}
GREY = '0.6'
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'axes.titlelocation': 'left'})
short = lambda g: g.replace('Pluripotent/Pluripotent-derived', 'Pluripotent').replace('Fetal life support', 'Fetal LS')
DEPTH = dict(cmap='RdBu_r', norm=mcolors.TwoSlopeNorm(0, -4, 2))


# ---------- shared pieces ----------

def sample_order(d):
    """Samples grouped in fit order, highest SPOT first within a group; returns order and group block bounds."""
    order, bounds = [], []
    for i in d.idx:
        bounds.append((len(order), len(order) + len(i)))
        order += list(i[np.argsort(-np.nan_to_num(d.spot[i]))])
    return np.array(order), bounds


def box(ax, s, e, y0, y1, kind, called, label=None):
    ls = '-' if called else (0, (2, 1.5))
    ax.add_patch(Rectangle((s - .5, y0), e - s, y1 - y0, fill=False, ec='white', lw=3.4 if called else 2.6, zorder=5))
    ax.add_patch(Rectangle((s - .5, y0), e - s, y1 - y0, fill=False, ec=CALL[kind], lw=1.8 if called else 1.3, ls=ls, zorder=6))
    if label:
        ax.text(e - .5, y0, label, fontsize=7, weight='bold', color=CALL[kind], ha='left', va='bottom', zorder=7,
                bbox=dict(fc='white', ec='none', pad=.4, alpha=.85))


def fp_lines(ax, d, **kw):
    for s, e in runs(d.fit.footprinted):
        for x in (s - .5, e - .5): ax.axvline(x, c='k', lw=kw.get('lw', .6), zorder=4)


def heat(ax, X, extent_rows, label, cbar_ticks=None, **kw):
    n = X.shape[1]
    im = ax.imshow(X, aspect='auto', interpolation='nearest', extent=(-.5, n - .5, extent_rows, 0), **kw)
    cb = plt.colorbar(im, cax=ax.inset_axes([1.012, .08, .008, .84]), ticks=cbar_ticks)
    cb.set_label(label, fontsize=8); cb.ax.tick_params(labelsize=7)
    for sp in ax.spines.values(): sp.set_visible(True)
    return im


def group_rows(ax, d):
    G = len(d.fit.group_names)
    ax.set_yticks(np.arange(G) + .5, [short(g) for g in d.fit.group_names], fontsize=7.5)
    for y in range(1, G): ax.axhline(y, c='white', lw=.5)


def sample_rows(ax, bounds, d):
    ax.set_yticks([(a + b) / 2 for a, b in bounds], [short(g) for g in d.fit.group_names], fontsize=7.5)
    for a, _ in bounds[1:]: ax.axhline(a, c='k', lw=.6)


def strip(ax, values, groups, names, color=None, hollow=None, err=None, flag=None, labels=None, focal=None, s=12):
    """Per-sample values by group; focal group in `color`, others grey; optional error bars and flagged labels."""
    rng = np.random.default_rng(0)
    for k, g in enumerate(names):
        i = np.flatnonzero(groups == g); x = k + rng.uniform(-.28, .28, len(i))
        c = color if g == focal else GREY
        if err is not None: ax.errorbar(x, values[i], yerr=err[i], fmt='none', ecolor=c, lw=.5, alpha=.5, zorder=2, rasterized=True)
        filled = np.ones(len(i), bool) if hollow is None else ~hollow[i]
        ax.scatter(x[filled], values[i][filled], s=s, c=c, lw=0, zorder=3, rasterized=True)
        ax.scatter(x[~filled], values[i][~filled], s=s, facecolors='none', edgecolors=c, lw=.8, zorder=3, rasterized=True)
        ax.plot([k - .35, k + .35], [np.median(values[i])] * 2, c='k', lw=1.2, zorder=4)
        if flag is not None:
            for j in np.flatnonzero(flag[i]): ax.annotate(labels[i[j]], (x[j], values[i[j]]), fontsize=6, xytext=(3, 0), textcoords='offset points', va='center')
    ax.set_xticks(range(len(names)), [short(g) for g in names], rotation=45, ha='right', fontsize=7.5)
    ax.set_xlim(-.6, len(names) - .4)


# ---------- page 1: overview ----------

def overview(d, C, pdf):
    F = d.fit; n, G, S = F.p_any.size, len(F.group_names), len(F.sample_names)
    order, bounds = sample_order(d)
    rows = [('fp', .9), ('obs', 3.0), ('oe', 3.0), ('oe5', 3.0), ('m', 3.0), ('sd', 1.5), ('gmean', 1.9), ('more', 1.9), ('less', 1.9)]
    fig = plt.figure(figsize=(16, 27))
    gs = fig.add_gridspec(len(rows), 2, height_ratios=[h for _, h in rows], width_ratios=[1, .035], hspace=.3, wspace=.15,
                          left=.08, right=.9, top=.945, bottom=.03)
    ax = {k: fig.add_subplot(gs[i, 0]) for i, (k, _) in enumerate(rows)}
    for k in list(ax)[1:]: ax[k].sharex(ax['fp'])
    nc = int(C.called.sum()) if len(C) else 0
    fig.suptitle(f"{d.name}   {d.interval}   ·   {S} samples, {G} groups   ·   {len(runs(F.footprinted))} footprints   ·   "
                 f"{nc} calls ({RULES[d.rule]}), {len(C) - nc} near calls", x=.08, y=.985, ha='left', fontsize=13, weight='bold')
    warn = dhs_flags(d)
    fig.text(.08, .972, '\n'.join('▲ ' + w for w in warn) or 'QC: no warnings', fontsize=8.5, va='top',
             color='#b35c00' if warn else '#555', weight='bold')

    a = ax['fp']; a.plot(np.arange(n), F.p_any, c='k', lw=1, drawstyle='steps-mid')
    a.axhline(F.cfg.p_footprint, c=GREY, lw=.7, ls='--'); a.set_ylim(0, 1.05); a.set_ylabel('P(any group\ndepleted)', fontsize=8)
    for s_, e_ in runs(F.footprinted):
        a.plot([s_ - .5, e_ - .5], [1.12, 1.12], transform=a.get_xaxis_transform(), c='k', lw=5, solid_capstyle='butt', clip_on=False)
    a.set_title('Footprints (bars; edges as vertical lines throughout): runs ≥ 5 bp where some group mean is < −log2 1.5 with probability ≥ 0.95',
                fontsize=9, pad=14)

    rel = d.obs / np.maximum(d.obs.mean(1, keepdims=True), 1e-9)
    a = ax['obs']; heat(a, rel[order], S, 'cuts / sample mean', [0, 1, 2, 4], cmap='Greys', norm=mcolors.PowerNorm(.5, 0, 4)); sample_rows(a, bounds, d)
    a.set_title('1 · Data: observed cuts per base, relative to the sample\'s mean per base (rows: group, then SPOT3 high→low)', fontsize=9)
    side = fig.add_subplot(gs[1, 1], sharey=a); cov = np.log10(d.obs.sum(1) + 1)[order]
    side.barh(np.arange(S) + .5, cov - cov.min() * .9, left=cov.min() * .9, height=1, color='0.45')
    side.set_xlabel('log10 cuts', fontsize=7); side.tick_params(labelleft=False, labelsize=6)

    a = ax['oe']; heat(a, np.log2((d.obs + .5) / (d.exp + .5))[order], S, 'log2 O/E per base', [-4, -2, 0, 2], **DEPTH); sample_rows(a, bounds, d)
    a.set_title('2 · Data over expected: log2((obs + ½) / (exp + ½)) per base, no smoothing', fontsize=9)
    side = fig.add_subplot(gs[2, 1], sharey=a)
    side.barh(np.arange(S) + .5, d.spot[order], height=1, color='0.45'); side.set_xlabel('SPOT3', fontsize=7)
    side.tick_params(labelleft=False, labelsize=6); side.set_xlim(0, max(.85, np.nanmax(d.spot) * 1.05))

    lr = lambda X: np.convolve(X, np.ones(5), 'same')
    oe5 = np.log2((np.apply_along_axis(lr, 1, d.obs) + .5) / (np.apply_along_axis(lr, 1, d.exp) + .5))
    a = ax['oe5']; heat(a, oe5[order], S, 'log2 O/E, 5-bp sums', [-4, -2, 0, 2], **DEPTH); sample_rows(a, bounds, d)
    a.set_title('3 · Same, smoothed: obs and exp summed over a centred 5-bp window before the ratio', fontsize=9)
    a = ax['m']; heat(a, F.theta_mean[order], S, 'segmented depth m (log2)', [-4, -2, 0, 2], **DEPTH); sample_rows(a, bounds, d)
    a.set_title('4 · Stage 1: each sample segmented on its own (posterior mean depth)', fontsize=9)
    a = ax['sd']; heat(a, np.sqrt(F.theta_var[order]), S, 'posterior SD', [0, .5, 1, 1.5], cmap='Purples', vmin=0, vmax=1.5); sample_rows(a, bounds, d)
    a.set_title('5 · Stage 1 uncertainty (posterior SD of the depth)', fontsize=9)
    a = ax['gmean']; heat(a, F.mean, G, 'group mean depth (log2)', [-4, -2, 0, 2], **DEPTH); group_rows(a, d)
    a.set_title(f'6 · Group mean depth; solid boxes = calls ({RULES[d.rule]}), dashed = best 5-bp window of a near call', fontsize=9)

    for kind in ('more', 'less'):
        a = ax[kind]
        cmap = mcolors.LinearSegmentedColormap.from_list(kind, ['white', mcolors.to_rgba(CALL[kind], .35), CALL[kind]])
        X = np.full((G, n), np.nan); w = d.wt['w']; X[:, w // 2:w // 2 + d.wt['delta'].shape[1]] = d.wclose(kind)
        heat(a, np.clip(X, 0, 1.3), G, 'closeness', [0, .5, .85, 1], cmap=cmap, vmin=0, vmax=1.3); group_rows(a, d)
        a.set_title(f'{7 if kind == "more" else 8} · Closeness of the 5-bp window centred on each base to a "{kind} depleted" call: '
                    f'min(lower bound of A / {F.cfg.a_min}, {"−" if kind == "more" else ""}Δ / log2 1.5), depths averaged over the window; ≥ 1 = window passes', fontsize=9)

    for _, c in C.iterrows():
        g = c.g; b0, b1 = bounds[g]
        spans = runs(d.called(c.kind)[g, c.start:c.end]) if c.called else [(c.win.min() - c.start, c.win.max() + 1 - c.start)]
        for i, (s_, e_) in enumerate(spans):
            s_, e_ = s_ + c.start, e_ + c.start
            for k in ('gmean', c.kind): box(ax[k], s_, e_, g, g + 1, c.kind, c.called, c.tag if i == 0 else None)
            for k in ('oe5', 'm'): box(ax[k], s_, e_, b0, b1, c.kind, c.called)
    for k in ax: fp_lines(ax[k], d, lw=.5)
    for k in list(ax)[:-1]: ax[k].tick_params(labelbottom=False)
    ax['less'].set_xlabel('position in DHS (bp)'); ax['fp'].set_xlim(-.5, n - .5)
    save(pdf, fig)


# ---------- candidate rows ----------

def span(bases):
    return f'{bases.min()}–{bases.max() + 1}'


def verdict_lines(c, cfg, d):
    ok = lambda b: '✓' if b else '✗'
    small, grey = dict(fontsize=7.5), dict(fontsize=7.5, color='#555')
    L = c.end - c.start
    lines = [(f"{c.tag}  {short(c.group)} · {c.kind} · {c.start}–{c.end}", dict(weight='bold', fontsize=10)),
             (('CALLED' if c.called else 'NOT CALLED: ' + c.outcome) + f'   closeness {c.close:.2f}', dict(weight='bold', color=CALL[c.kind] if c.called else '#555')),
             (f"best 5-bp window {span(c.win)}, depths averaged over it:", grey),
             (f"{ok(c.lo >= cfg.a_min)} A {c.a:.2f}, lower bound {c.lo:.2f}; needs ≥ {cfg.a_min}", small),
             (f"{ok(SIGN[c.kind] * c.delta > cfg.fold)} Δ {c.delta:+.2f}; needs {'<' if c.kind == 'more' else '>'} {SIGN[c.kind] * cfg.fold:+.3f}", small),
             (f"bases passing on their own: {c.n_pass_win}/5 in the window, {c.n_pass}/{L} in the run", dict(fontsize=7.5, weight='bold')),
             (f"called: {c.n_called} bp ({RULES[d.rule]});  {RULES[d.alt]}: {c.alt_bp} bp", small),
             ('each base of the window on its own:', grey)] + [('   ' + t, dict(fontsize=7, color='#333')) for t in c.detail]
    lines += [(f"ties (π < {cfg.a_min}): {short(', '.join(c.ties.split(', ')[:2]))};  side {c.n_side}/{c.n_g}", small),
              (f"expected cuts/bp: group {c.exp_g:.0f}, others {c.exp_o:.0f};  posterior SD {c.sd_g:.2f}, between-sample SD {c.spread_g:.2f}", small),
              (f"scaling slope {c.slope:.2f} predicts Δ {c.delta_scale:+.2f}", small)]
    return lines + [('▲ ' + f, dict(fontsize=7.5, color='#b35c00', weight='bold')) for f in c.warn]


def candidate_row(fig, sgs, d, c):
    F, cfg = d.fit, d.fit.cfg; n = F.p_any.size; g, kind = c.g, c.kind; col = CALL[kind]; win = c.win
    inner = sgs.subgridspec(2, 5, width_ratios=[3.0, 3.2, 1.7, 1.8, 3.3], height_ratios=[2.2, 1], wspace=.42, hspace=.12)
    lo_, hi_ = max(0, c.start - 15), min(n, c.end + 15); x = np.arange(lo_, hi_)
    others = d.others[g]

    a = fig.add_subplot(inner[0, 0])
    for h in range(len(F.group_names)):
        if h != g: a.plot(x, F.mean[h, x], c='0.8', lw=.7, drawstyle='steps-mid')
    a.fill_between(x, F.mean[g, x] - F.se[g, x], F.mean[g, x] + F.se[g, x], step='mid', color=col, alpha=.25, lw=0)
    a.plot(x, F.mean[g, x], c=col, lw=1.8, drawstyle='steps-mid')
    a.plot(x, others[x], c='k', lw=1, ls='--', drawstyle='steps-mid')
    a.plot(x, others[x] + SIGN[kind] * cfg.fold, c='k', lw=.6, ls=':', drawstyle='steps-mid')
    a.axvspan(win.min() - .5, win.max() + .5, color='k', alpha=.08, lw=0)
    top = lambda s_, e_, cc, y, lw: a.plot([s_ - .5, e_ - .5], [y, y], transform=a.get_xaxis_transform(), c=cc, lw=lw, solid_capstyle='butt', clip_on=False)
    for s_, e_ in runs(d.called(kind)[g, lo_:hi_]): top(lo_ + s_, lo_ + e_, col, 1.02, 4)
    for s_, e_ in runs(base_pass(F, kind)[g, lo_:hi_]): top(lo_ + s_, lo_ + e_, '0.3', 1.075, 1.6)
    for s_, e_ in runs(F.footprinted):
        for xx in (s_ - .5, e_ - .5):
            if lo_ <= xx <= hi_: a.axvline(xx, c='k', lw=.5)
    a.set_ylabel('group mean\ndepth (log2)', fontsize=8); a.tick_params(labelbottom=False); a.set_xlim(lo_ - .5, hi_ - .5)
    a.set_title(c.tag, fontsize=11, weight='bold', color=col, pad=8)

    b = fig.add_subplot(inner[1, 0], sharex=a)
    A, LB = getattr(F, f'a_{kind}')[g, x], getattr(F, f'lo_{kind}')[g, x]
    b.fill_between(x, LB, A, step='mid', color=col, alpha=.25, lw=0); b.plot(x, A, c=col, lw=1.2, drawstyle='steps-mid')
    b.plot(x, LB, c=col, lw=.6, drawstyle='steps-mid')
    w = d.wt['w']; t = np.arange(max(0, lo_ - w // 2), min(d.wt['delta'].shape[1], hi_ - w // 2))
    b.plot(t + w // 2, d.wt[f'lo_{kind}'][g, t], c='k', lw=1, drawstyle='steps-mid')
    b.axhline(cfg.a_min, c='k', lw=.7, ls='--'); b.set_ylim(.4, 1); b.axvspan(win.min() - .5, win.max() + .5, color='k', alpha=.08, lw=0)
    b.set_ylabel('A (to lower\nbound)', fontsize=8); b.set_xlabel('bp', fontsize=8)

    a = fig.add_subplot(inner[:, 1])
    mb, sd = F.theta_mean[:, win].mean(1), np.sqrt(F.theta_var[:, win].mean(1))
    lab = np.array(F.sample_groups); q = np.full(len(mb), 1.); q[d.idx[g]] = sample_side(d, g, kind, win)
    strip(a, mb, lab, F.group_names, color=col, hollow=q <= .5, err=sd, focal=c.group, s=10)
    o = np.mean(others[win])
    a.axhline(o, c='k', lw=.8, ls='--'); a.axhspan(o, o + SIGN[kind] * cfg.fold, color=col, alpha=.07, lw=0)
    a.plot([g - .45, g + .45], [o + c.delta_scale] * 2, c=col, lw=1.5, ls=(0, (2, 1)))
    a.set_ylabel(f'sample depth over {span(win)}\n(± posterior SD)', fontsize=8)
    a.set_title('samples (hollow = not on the called side)', fontsize=7.5)

    a = fig.add_subplot(inner[:, 2])
    pi = pairwise(d, g, kind, win)
    a.barh(range(len(pi)), pi.values, color=[GREY if v >= cfg.a_min else '#c9a227' for v in pi.values], height=.7)
    a.barh(len(pi) + .4, c.a, color=col, height=.9)
    a.errorbar(c.lo, len(pi) + .4, xerr=[[0], [c.a - c.lo]], fmt='|', c='k', ms=8, lw=1)
    a.set_yticks(list(range(len(pi))) + [len(pi) + .4], [short(h) for h in pi.index] + ['A (lower bound)'], fontsize=7)
    a.axvline(cfg.a_min, c='k', lw=.7, ls='--'); a.axvline(.5, c=GREY, lw=.5); a.set_xlim(0, 1)
    a.set_xlabel('π: P(on the called side)', fontsize=8); a.set_title('vs each group', fontsize=8)

    a = fig.add_subplot(inner[:, 3]); other = lab != c.group
    a.scatter(d.spot[other], mb[other], s=7, c=GREY, lw=0, rasterized=True); a.scatter(d.spot[~other], mb[~other], s=12, c=col, lw=0, rasterized=True)
    a.axhline(o, c='k', lw=.8, ls='--'); a.set_xlabel('SPOT3', fontsize=8); a.set_title('depth vs sample quality', fontsize=8)

    a = fig.add_subplot(inner[:, 4]); a.axis('off'); y = 1.02
    for text, kw in verdict_lines(c, cfg, d):
        kw = dict(kw); fs = kw.pop('fontsize', 8.5)
        for k, part in enumerate(textwrap.wrap(text, int(64 * 7.5 / fs), subsequent_indent='     ')):
            a.text(-.08, y, part, transform=a.transAxes, va='top', fontsize=fs, **kw); y -= .062 * fs / 7.5


def legend(rule):
    calls = ('Calls (per-base rule): bases where the lower bound of A ≥ 0.75 and |Δ| > log2 1.5, in runs ≥ 5 bp.  '
             if rule == 'base' else
             'Calls (window rule): depths are averaged over each 5-bp window and the window is tested (lower bound of A ≥ 0.75 and |Δ| > log2 1.5); the best passing '
             'window is called and extended through flanking bases that pass on their own, then the next best window not touching a call, and so on.  ')
    return (calls + 'Each candidate is reported on its best 5-bp window (grey band), with depths averaged over it.\n'
            'Left: group mean (colour) ± SE, other groups (light grey), median of the other groups (dashed), 1.5-fold from it (dotted); coloured bar on top = called bases, '
            'thin dark bar = bases passing on their own; black verticals = footprint edges.  Below it: A per base shaded down to its lower bound (colour), '
            'lower bound of the 5-bp window centred there (black); dashed = 0.75.\n'
            'Middle: each sample\'s depth averaged over the window ± posterior SD; black ticks = group medians; shaded = from the others\' median to 1.5-fold beyond it; '
            'dashed coloured tick = the depth group-wide scaling alone predicts.  Bars: π against each other group on the window (gold = tie, < 0.75), and A with its lower bound.  '
            'Right: the same sample depths against SPOT3.')


def candidate_pages(d, C, pdf, per_page=6):
    for p in range(0, len(C), per_page):
        chunk = C.iloc[p:p + per_page]
        H = 3.0 * per_page + 2.1
        fig = plt.figure(figsize=(16, H))
        gs = fig.add_gridspec(per_page, 1, hspace=.45, left=.06, right=.99, top=1 - 2.05 / H, bottom=.03)
        fig.suptitle(f'{d.name} · candidates {p + 1}–{p + len(chunk)} of {len(C)}', x=.06, y=1 - .25 / H, ha='left', fontsize=12, weight='bold')
        fig.text(.06, 1 - .5 / H, '\n'.join(textwrap.fill(p, 215) for p in legend(d.rule).split('\n')), fontsize=7.5, va='top', color='#333', linespacing=1.4)
        for i, (_, c) in enumerate(chunk.iterrows()): candidate_row(fig, gs[i], d, c)
        save(pdf, fig)


# ---------- QC page ----------

def qc_page(d, pdf):
    F = d.fit; qc, pit = d.qc; names = list(F.group_names); lab = qc.group.to_numpy(); labels = qc['sample'].to_numpy()
    fig, axs = plt.subplots(2, 4, figsize=(16, 9), gridspec_kw=dict(hspace=.62, wspace=.32, left=.05, right=.98, top=.9, bottom=.14))
    fig.suptitle(f'{d.name} · checks of the assumptions', x=.05, ha='left', fontsize=12, weight='bold')

    a = axs[0, 0]; bins = np.linspace(0, 1, 21)
    a.hist(pit['pit_open'], bins=bins, density=True, color='0.75', label='outside footprints')
    a.hist(pit['pit_fp'], bins=bins, density=True, histtype='step', color='k', lw=1.2, label='in footprints')
    a.axhline(1, c=CALL['less'], lw=.8, ls='--'); a.legend(fontsize=7, frameon=False, loc='upper center')
    a.set_xlabel('randomized PIT of each count'); a.set_ylabel('density')
    a.set_title('Count model calibration\n(flat = NB with the sample\'s dispersion fits)', fontsize=9)

    a = axs[0, 1]; strip(a, qc.disp.to_numpy(), lab, names, flag=qc.disp.to_numpy() > 2.5, labels=labels)
    a.axhline(1, c=CALL['less'], lw=.8, ls='--'); a.set_ylabel('mean squared Pearson residual')
    a.set_title('Per-sample dispersion index\n(1 = as noisy as the model expects)', fontsize=9)

    a = axs[0, 2]; off = qc.offset.to_numpy(); strip(a, off, lab, names, flag=np.abs(off - np.median(off)) > 1, labels=labels)
    a.axhline(0, c=CALL['less'], lw=.8, ls='--'); a.set_ylabel('median depth outside footprints')
    a.set_title('Per-sample offset\n(open DNA assumed at 0)', fontsize=9)

    a = axs[0, 3]; loo = qc.loo.to_numpy()
    med = pd.Series(loo).groupby(lab).transform('median').to_numpy()
    worst = loo <= np.sort(loo)[min(5, len(loo) - 1)]
    strip(a, loo, lab, names, flag=(loo < med - .2) & worst, labels=labels); a.set_ylabel('correlation with own group')
    a.set_title('Within-group consistency\n(profile vs rest of its group)', fontsize=9)

    a = axs[1, 0]; a.scatter(qc.spot, qc.fp_depth, s=8, c=GREY, lw=0)
    for gname, sub in qc.groupby('group'):
        a.annotate(short(gname), (sub.spot.median(), sub.fp_depth.median()), fontsize=6.5, ha='center', weight='bold')
    ok = qc[['spot', 'fp_depth']].dropna(); r = np.corrcoef(ok.spot, ok.fp_depth)[0, 1]
    a.set_xlabel('SPOT3'); a.set_ylabel('mean depth in footprints')
    a.set_title(f'Footprint depth vs sample quality\n(r = {r:.2f}; labels at group medians)', fontsize=9)

    a = axs[1, 1]; excl = candidate_mask(d); sl = np.array([scaling(d, g, excl[g])[0] for g in range(len(names))])
    a.bar(range(len(names)), sl, color=['0.25' if abs(v - 1) > .15 else '0.7' for v in sl], width=.7)
    a.axhline(1, c='k', lw=.7); a.set_xticks(range(len(names)), [short(g) for g in names], rotation=45, ha='right', fontsize=7.5)
    a.set_ylabel('slope vs median of others'); a.set_title('Group-wide depth scaling in footprints\n(robust slope, own candidates excluded; dark = off by > 15%)', fontsize=9)

    a = axs[1, 2]; strip(a, np.log10(qc.cuts.to_numpy() + 1), lab, names)
    a.set_ylabel('log10 cuts in DHS'); a.set_title('Coverage per sample', fontsize=9)

    a = axs[1, 3]; fp = F.footprinted
    for k, i in enumerate(d.idx):
        spread = np.median(F.theta_mean[i][:, fp].std(0)); noise = np.median(np.sqrt(F.theta_var[i][:, fp]))
        a.scatter(noise, spread, s=14, c='0.3'); a.annotate(short(names[k]), (noise, spread), fontsize=6.5, xytext=(3, 2), textcoords='offset points')
    lim = a.get_xlim()[1], a.get_ylim()[1]; m_ = max(lim); a.plot([0, m_], [0, m_], c=GREY, lw=.7, ls='--')
    a.set_xlim(0, m_); a.set_ylim(0, m_)
    a.set_xlabel('count noise: median posterior SD'); a.set_ylabel('between-sample SD')
    a.set_title('Heterogeneity vs count noise in footprints\n(above diagonal = samples differ beyond noise)', fontsize=9)
    save(pdf, fig)


def table_page(d, C, pdf):
    cols = ['tag', 'group', 'kind', 'run', 'best window', 'outcome', 'closeness', 'A', 'lower bd', 'Δ', 'pass in window', 'pass in run',
            'called bp', f'{RULES[d.alt]} bp', 'ties (π < 0.75)', 'side', 'warnings']
    rows = [[c.tag, short(c.group), c.kind, f'{c.start}–{c.end}', span(c.win), 'called' if c.called else c.outcome, f'{c.close:.2f}',
             f'{c.a:.2f}', f'{c.lo:.2f}', f'{c.delta:+.2f}', f'{c.n_pass_win}/5', f'{c.n_pass}/{c.end - c.start}', c.n_called or '', c.alt_bp or '',
             short(c.ties), f'{c.n_side}/{c.n_g}', '; '.join(c.warn)] for _, c in C.iterrows()]
    H = 1.2 + .24 * max(len(rows), 1)
    fig = plt.figure(figsize=(16, H)); a = fig.add_axes([.02, .02, .96, 1 - .9 / H]); a.axis('off')
    fig.suptitle(f'{d.name} · candidates ({RULES[d.rule]}; called first, then near calls by closeness); A, lower bound and Δ are for the best 5-bp window',
                 x=.02, y=1 - .2 / H, ha='left', fontsize=12, weight='bold')
    if rows:
        t = a.table(cellText=rows, colLabels=cols, loc='upper left', cellLoc='left', colLoc='left',
                    colWidths=[.025, .06, .03, .045, .05, .11, .045, .03, .045, .035, .055, .05, .04, .06, .16, .03, .135])
        t.auto_set_font_size(False); t.set_fontsize(7.5); t.scale(1, 1.25)
        for (i, j), cell in t.get_celld().items():
            cell.set_edgecolor('0.85'); cell.set_linewidth(.4)
            if i == 0: cell.set_text_props(weight='bold')
            elif j == 2: cell.get_text().set_color(CALL[rows[i - 1][2]])
            elif j == 5 and rows[i - 1][5] == 'called': cell.get_text().set_weight('bold')
            elif j == 16: cell.get_text().set_color('#b35c00')
    save(pdf, fig)


def save(pdf, fig):
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def report(d, path, near=.85):
    """Write the report PDF; returns the candidate table."""
    C = candidates(d, near=near)
    with PdfPages(path) as pdf:
        overview(d, C, pdf); table_page(d, C, pdf); candidate_pages(d, C, pdf); qc_page(d, pdf)
    return C
