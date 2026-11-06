"""Window rule: test pooled 5-bp windows; call the best passing window and extend it through bases that pass on their own."""
import numpy as np

from footprint_tools.stats.differential_bayesian.calls import discriminability, group_tracks

SIGN = {'more': -1, 'less': 1}


def window_tracks(F, idx, w=5):
    """A, lower bound and Delta of every w-bp window (index t = first base) on window-averaged depths.
    Variances are averaged too, i.e. bases are treated as fully correlated (conservative)."""
    k = np.ones(w) / w
    avg = lambda X: np.apply_along_axis(lambda r: np.convolve(r, k, 'valid'), 1, X.astype(float))
    m, v = avg(F.theta_mean), avg(F.theta_var)
    _, _, delta = group_tracks(m, v, idx)
    return dict(delta=delta, w=w, **discriminability(m, v, idx, F.cfg.n_boot, F.cfg.seed))


def base_pass(F, kind):
    c = F.cfg
    return (getattr(F, f'lo_{kind}') >= c.a_min) & (SIGN[kind] * F.delta > c.fold)


def window_pass(F, T, kind):
    c = F.cfg
    return (T[f'lo_{kind}'] >= c.a_min) & (SIGN[kind] * T['delta'] > c.fold)


def window_closeness(F, T, kind):
    c = F.cfg
    return np.minimum(T[f'lo_{kind}'] / c.a_min, SIGN[kind] * T['delta'] / c.fold)


def cover(mask_w, n, w):
    """Bases covered by the windows set in mask_w (rows x windows) -> rows x n."""
    out = np.zeros((mask_w.shape[0], n), bool)
    for i in range(w): out[:, i:i + mask_w.shape[1]] |= mask_w
    return out


def called_window(F, T, kind):
    """G x n. Greedy: take the best passing window not touching a call, call it and extend it through contiguous
    flanking bases that pass on their own; repeat until no passing window is left."""
    ok, close, base = window_pass(F, T, kind), window_closeness(F, T, kind), base_pass(F, kind)
    G, n, w = ok.shape[0], F.p_any.size, T['w']
    out = np.zeros((G, n), bool)
    for g in range(G):
        for t in sorted(np.flatnonzero(ok[g]), key=lambda t: -close[g, t]):
            if out[g, t:t + w].any(): continue
            s, e = t, t + w
            while s > 0 and base[g, s - 1] and not out[g, s - 1]: s -= 1
            while e < n and base[g, e] and not out[g, e]: e += 1
            out[g, s:e] = True
    return out
