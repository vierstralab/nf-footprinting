
"""Per-base and per-DHS p-values for differential footprint occupancy, from the stage-1 per-sample posteriors.

Effect: W_S, the probability-weighted share of the expected depth difference between the samples of a set S of groups
and the rest that lies on S's side (pairs weighted by E|m_t - m_s|). Statistic per base: the largest |W_S - 1/2| - e,
over every set of up to 5 groups (both directions), divided by its grouped delete-one jackknife SD. Reference: the same
statistic under label permutations (studentized permutation max-T); one permutation is shared by all bases, so the DHS
statistics keep the dependence between bases. Beyond the permutations, a generalised Pareto tail fitted to the top
values (Knijnenburg et al. 2009). e = 0 tests "no difference" (exact under exchangeable labels); e = 0.1, 0.2, 0.3 test
W below 0.6, 0.7, 0.8 and reuse the same reference for the studentized shortfall, so they are approximate.
"""
import argparse, itertools

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.special import ndtr

K_MAX, N_EXC_FRAC, MIN_HITS, WIN, BATCH = 5, .1, 10, 5, 40
EFFECTS = {'any': 0., 'w60': .1, 'w70': .2, 'w80': .3}       # null: |W - 1/2| <= e


def pair_terms(m, v):
    """pos[j, s, t] = E[(m_t - m_s)_+] (s more depleted than t) and ab[j, s, t] = E|m_t - m_s|, positions x N x N."""
    mu = m.T[:, None, :] - m.T[:, :, None]; sd = np.sqrt(v.T[:, :, None] + v.T[:, None, :]); z = mu / sd
    pos = mu * ndtr(z) + sd * np.exp(-z ** 2 / 2) / np.sqrt(2 * np.pi)
    return pos.astype(np.float32), (2 * pos - mu).astype(np.float32)


def set_weights(G):
    """All sets of up to K_MAX groups as weights a_S with W_S = a_S . R - off_S, R_g = sum over h of W_gh."""
    S = [s for k in range(1, K_MAX + 1) for s in itertools.combinations(range(G), k)]
    ind = np.zeros((len(S), G))
    for i, s in enumerate(S): ind[i, list(s)] = 1
    k = ind.sum(1)
    return S, ind / (k * (G - k))[:, None], (k - 1) / (2 * (G - k))


def group_sums(pos, ab, labelings, G):
    """For each labeling: sums over t in group h of pos[s, t] and ab[s, t] (N x G x positions), in one wide product."""
    n, N = pos.shape[:2]
    Z = np.concatenate([np.eye(G, dtype=np.float32)[g] for g in labelings], 1)
    out = [(x.reshape(-1, N) @ Z).reshape(n, N, len(labelings), G).transpose(2, 1, 3, 0) for x in (pos, ab)]
    return zip(labelings, *out)


def set_stats(g, PZ, AZ, G, a, off):
    """W_S and its grouped delete-one jackknife SD for every set at every position (sets x positions)."""
    Z = np.eye(G)[g]; ns = Z.sum(0); N, n = len(g), PZ.shape[2]
    num, den = [(Z.T @ X.reshape(N, -1)).reshape(G, G, n) for X in (PZ, AZ)]
    off_diag = ~np.eye(G, dtype=bool)[..., None]
    W = np.where(off_diag, num / np.where(off_diag, den, 1), 0)
    dW = (num[g] - PZ) / (den[g] - AZ) - W[g]                      # change in W_{g(s) h} when sample s is deleted
    dW[np.arange(N), g] = 0
    D = -dW; D[np.arange(N), g] = dW.sum(1)
    for h in range(G):
        i = g == h; D[i] = np.sqrt((ns[h] - 1) / ns[h]) * (D[i] - D[i].mean(0))
    C = np.einsum('ngj,nhj->ghj', D, D)
    var = ((a @ C.reshape(G, -1)).reshape(len(a), G, n) * a[:, :, None]).sum(1)
    return a @ W.sum(1) - off[:, None], np.sqrt(np.maximum(var, 1e-12))


def log10_p(Mrep, T):
    """log10 P(M >= T) per column: (1 + hits) / (B + 1) when at least MIN_HITS permutations reach T, else a generalised
    Pareto tail fitted to the top 10% by probability-weighted moments (Hosking & Wallis 1987), xi < 0 used as 0.
    Also returns the fitted xi."""
    B = len(Mrep); n_exc = max(int(N_EXC_FRAC * B), 20)
    s = np.sort(Mrep, 0); u = s[-n_exc - 1].astype(float); y = s[-n_exc:] - u
    i = np.arange(n_exc)[:, None]
    a0, a1 = y.mean(0), (y * (n_exc - 1 - i) / (n_exc - 1)).mean(0)
    xi_fit = 2 - a0 / np.maximum(a0 - 2 * a1, 1e-12)
    xi = np.clip(xi_fit, 0, None); scale = np.maximum(np.where(xi > 0, 2 * a0 * a1 / np.maximum(a0 - 2 * a1, 1e-12), a0), 1e-9)
    x = np.maximum(T - u, 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        surv = np.where(xi > 0, -np.log1p(xi * x / scale) / np.where(xi > 0, xi, 1), -x / scale) / np.log(10)
    hits = (Mrep >= T).sum(0)
    return np.where(hits >= MIN_HITS, np.log10((1 + hits) / (B + 1)), np.log10(n_exc / B) + surv), xi_fit


def pvalues(m, v, labels, n_perm, rng):
    names = list(dict.fromkeys(labels)); g = np.array([names.index(x) for x in labels]); G = len(names)
    _, a, off = set_weights(G); pos, ab = pair_terms(m, v)
    (_, PZ, AZ), = group_sums(pos, ab, [g], G); W, SD = set_stats(g, PZ, AZ, G, a, off)
    M = []
    for b in range(0, n_perm, BATCH):
        for gp, PZ, AZ in group_sums(pos, ab, [rng.permutation(g) for _ in range(min(BATCH, n_perm - b))], G):
            w, sd = set_stats(gp, PZ, AZ, G, a, off); M.append((np.abs(w - .5) / sd).max(0))
    M = np.stack(M); Mw = sliding_window_view(M, WIN, axis=1).min(2).max(1, keepdims=True)
    out, cols = dict(names=np.array(names), c95=np.quantile(M, .95, 0).astype(np.float32)), np.arange(W.shape[1])
    for k, e in EFFECTS.items():
        t = (np.abs(W - .5) - e) / SD; i = t.argmax(0); T = t[i, cols]
        out[f'lp_{k}'], xi = (x.astype(np.float32) for x in log10_p(M, T))
        out[f'set_{k}'], out[f'w_{k}'] = i.astype(np.int16), W[i, cols].astype(np.float32)
        out[f'dhs_{k}'] = log10_p(M.max(1, keepdims=True), np.array([T.max()]))[0][0]
        out[f'dhs5_{k}'] = log10_p(Mw, np.array([sliding_window_view(T, WIN).min(1).max()]))[0][0]   # best 5-bp window
    out['xi'] = xi
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('dhs_ids', help='comma-separated'); p.add_argument('diff_data', nargs='+', help='per_dhs/diff_data.<dhs_id>.npz, same order')
    p.add_argument('--n_perm', type=int, default=1000); p.add_argument('--seed', type=int, default=0)
    p.add_argument('--summary', default='diff_pvalues.tsv')
    a = p.parse_args()
    rows = []
    for dhs_id, f in zip(a.dhs_ids.split(','), a.diff_data):
        with np.load(f) as z:
            m, v, lab, interval = z['theta_mean'].astype(float), z['theta_var'].astype(float), z['sample_groups'].astype(str), str(z['__interval__'])
        out = pvalues(m, v, lab, a.n_perm, np.random.default_rng(a.seed))
        np.savez_compressed(f'diff_pvalues.{dhs_id}.npz', interval=interval, **out)
        rows.append([dhs_id, interval, m.shape[1]] + [f"{out[f'{d}_{k}']:.3f}" for d in ('dhs', 'dhs5') for k in EFFECTS]
                    + [f"{out[f'lp_{k}'].min():.3f}" for k in EFFECTS])
    with open(a.summary, 'w') as fh:
        fh.write('\t'.join(['dhs_id', 'interval', 'width'] + [f'log10p_{d}_{k}' for d in ('dhs', 'dhs5') for k in EFFECTS]
                           + [f'min_log10p_base_{k}' for k in EFFECTS]) + '\n')
        fh.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
