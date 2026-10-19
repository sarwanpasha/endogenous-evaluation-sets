"""Controlled synthetic study of endogenous evaluation-set bias.

A rank-r ground truth F is generated, a reference prior L of controlled accuracy
is attached to it, and a one-parameter family of designs spends a fraction
alpha of its budget inside the model-defined band B and the rest uniformly over
the support.  alpha = 0 is the non-targeting control; alpha = 1 is the design
whose allocation coincides exactly with B.  Every design is scored on both the
model-free band B* and the model-defined band B, and we report the improvement
over the alpha = 0 control on each and the ratio of the two.

Nothing here is specific to drug screening: the object is a matrix, a
threshold, a prior, and a budget.
"""
import argparse, json, os, sys, time
import numpy as np

A = argparse.ArgumentParser()
A.add_argument('--P', type=int, default=800)          # rows
A.add_argument('--m', type=int, default=80)           # columns
A.add_argument('--rank', type=int, default=8)         # rank of the truth
A.add_argument('--eta', type=float, default=0.50)     # prior rms error, x sigma_F
A.add_argument('--sig', type=float, default=0.50)     # replicate sd, x sigma_F
A.add_argument('--q', type=float, default=0.90)       # threshold quantile
A.add_argument('--b', type=float, default=0.30)       # band half-width, x sigma_F
A.add_argument('--frac', type=float, default=0.32)    # budget as T/N
A.add_argument('--seeds', type=int, default=100)
A.add_argument('--alphas', type=str, default='0,0.25,0.5,0.75,1.0')
A.add_argument('--sweep', type=str, default='base')   # base|band|eta|rank|budget
A.add_argument('--out', type=str, default='syn_out.json')
A = A.parse_args()

ALPHAS = [float(x) for x in A.alphas.split(',')]


def make_truth(rs, P, m, r):
    """Rank-r matrix standardised to unit entry sd."""
    U = rs.standard_normal((P, r))
    V = rs.standard_normal((m, r))
    F = (U @ V.T) / np.sqrt(r)
    F = (F - F.mean()) / F.std()
    return F.ravel()


def run_one(rs, P, m, r, eta, sig, q, b, frac):
    """One replication.  Returns err[alpha] on B* and on B."""
    F = make_truth(rs, P, m, r)
    N = F.size
    tau = float(np.quantile(F, q))
    H = F >= tau

    # Reference prior: the truth seen through a model of controlled accuracy.
    # eta is its rms error in units of sigma_F (= 1 by construction).
    L = F + eta * rs.standard_normal(N)

    tstar = np.abs(F - tau) <= b          # model-free band  B*
    tband = np.abs(L - tau) <= b          # model-defined band B
    if tstar.sum() < 10 or tband.sum() < 10:
        return None

    T = max(int(frac * N), 1)
    lam = (sig ** 2) / (eta ** 2)         # exact, so estimation noise is not a confound
    idx_all = np.arange(N)
    idx_band = np.where(tband)[0]

    out = {}
    for al in ALPHAS:
        n_in = int(round(al * T))
        c = np.zeros(N, dtype=np.int64)
        s = np.zeros(N)
        for pool, n in ((idx_band, n_in), (idx_all, T - n_in)):
            if n <= 0 or pool.size == 0:
                continue
            pick = rs.integers(0, pool.size, size=n)
            cnt = np.bincount(pick, minlength=pool.size)
            hit = cnt > 0
            j = pool[hit]
            c[j] += cnt[hit]
            s[j] += cnt[hit] * F[j] + sig * np.sqrt(cnt[hit]) * rs.standard_normal(j.size)
        seen = c > 0
        Fh = L.copy()
        Fh[seen] = (lam * L[seen] + s[seen]) / (lam + c[seen])
        Hh = Fh >= tau
        out[al] = (float((Hh[tstar] != H[tstar]).mean()),
                   float((Hh[tband] != H[tband]).mean()))
    return out, int(tstar.sum()), int(tband.sum()), float((tstar & tband).mean() * N)


def sweep_settings():
    if A.sweep == 'band':
        return [('b', v, dict(b=v)) for v in (0.05, 0.10, 0.20, 0.30, 0.50)]
    if A.sweep == 'eta':
        return [('eta', v, dict(eta=v)) for v in (0.125, 0.25, 0.50, 1.00, 2.00)]
    if A.sweep == 'rank':
        return [('rank', v, dict(rank=v)) for v in (2, 4, 8, 16, 32)]
    if A.sweep == 'budget':
        return [('frac', v, dict(frac=v)) for v in (0.04, 0.16, 0.32, 0.64, 1.28)]
    return [('base', 0.0, {})]


res = []
t0 = time.time()
for name, val, over in sweep_settings():
    kw = dict(P=A.P, m=A.m, r=A.rank, eta=A.eta, sig=A.sig, q=A.q, b=A.b, frac=A.frac)
    if 'rank' in over:
        kw['r'] = over['rank']
    for k in ('b', 'eta', 'frac'):
        if k in over:
            kw[k] = over[k]
    acc = {al: [[], []] for al in ALPHAS}
    ms, mb, ov = [], [], []
    for sd in range(A.seeds):
        rs = np.random.default_rng(7717 * (1 + int(val * 1000)) + 31 * sd + hash(name) % 997)
        got = run_one(rs, **kw)
        if got is None:
            continue
        o, a, bb, c = got
        ms.append(a); mb.append(bb); ov.append(c)
        for al in ALPHAS:
            acc[al][0].append(o[al][0])
            acc[al][1].append(o[al][1])
    if not ms:
        continue
    row = {'sweep': name, 'value': val, 'setting': {k: (v if not isinstance(v, np.ndarray) else None)
                                                    for k, v in kw.items()},
           'Mstar': float(np.mean(ms)), 'Mband': float(np.mean(mb)),
           'overlap': float(np.mean(ov)), 'seeds': len(ms), 'arms': {}}
    base_s = float(np.mean(acc[0.0][0])) if 0.0 in acc else None
    base_b = float(np.mean(acc[0.0][1])) if 0.0 in acc else None
    for al in ALPHAS:
        es = np.array(acc[al][0]); eb = np.array(acc[al][1])
        gs = base_s - es.mean() if base_s is not None else None
        gb = base_b - eb.mean() if base_b is not None else None
        row['arms'][str(al)] = {
            'Bstar': round(float(es.mean()), 5), 'B': round(float(eb.mean()), 5),
            'Bstar_se': round(float(es.std(ddof=1) / np.sqrt(es.size)), 5),
            'B_se': round(float(eb.std(ddof=1) / np.sqrt(eb.size)), 5),
            'gain_Bstar': None if gs is None else round(float(gs), 5),
            'gain_B': None if gb is None else round(float(gb), 5),
            'x': (round(float(gb / gs), 4) if gs is not None and gs > 2e-3 and gb > 0 else None)}
    res.append(row)
    print(json.dumps({'k': 'tick', 'sweep': name, 'value': val,
                      'x@1.0': row['arms'][str(ALPHAS[-1])]['x'],
                      'sec': round(time.time() - t0, 1)}), flush=True)

json.dump(res, open(A.out + '.tmp', 'w'), indent=1)
os.replace(A.out + '.tmp', A.out)
print(json.dumps({'k': 'done', 'rows': len(res), 'out': A.out,
                  'sec': round(time.time() - t0, 1)}), flush=True)
