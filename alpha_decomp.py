"""Decomposition of the inflation factor in the controlled model.

For every targeting strength alpha the script measures the three quantities the
ratio proposition names,

    omega = |B n B*| / |B*|,
    kappa = mean gain on B n B* / mean gain on B,
    mu    = mean gain on B* \\ B / mean gain on B,

where a gain is the drop in boundary error relative to the alpha = 0 control on
the stated set, and compares the predicted factor 1 / (omega kappa + (1 - omega)
mu) with the factor measured directly.  The point of the exercise is that the
factor need not move monotonically in alpha: targeting raises the gain on the
overlap, which raises kappa, and at the same time leaves the part of the
model-free band the model misses untouched, which drives mu toward zero, and the
two act on the denominator in opposite directions.

Same generative model, defaults and seed convention as syn.py, so the alpha = 0
and alpha = 1 arms here reproduce the numbers reported there.
"""
import argparse, json, os, time
import numpy as np

A = argparse.ArgumentParser()
A.add_argument('--P', type=int, default=800)
A.add_argument('--m', type=int, default=80)
A.add_argument('--rank', type=int, default=8)
A.add_argument('--eta', type=float, default=0.50)
A.add_argument('--sig', type=float, default=0.50)
A.add_argument('--q', type=float, default=0.90)
A.add_argument('--b', type=float, default=0.30)
A.add_argument('--frac', type=float, default=0.32)
A.add_argument('--seeds', type=int, default=120)
A.add_argument('--alphas', type=str, default='0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0')
A.add_argument('--out', type=str, default='alpha_decomp.json')
A = A.parse_args()

ALPHAS = [float(x) for x in A.alphas.split(',')]


def make_truth(rs, P, m, r):
    U = rs.standard_normal((P, r))
    V = rs.standard_normal((m, r))
    F = (U @ V.T) / np.sqrt(r)
    return ((F - F.mean()) / F.std()).ravel()


def run_one(rs):
    F = make_truth(rs, A.P, A.m, A.rank)
    N = F.size
    tau = float(np.quantile(F, A.q))
    H = F >= tau
    L = F + A.eta * rs.standard_normal(N)
    tstar = np.abs(F - tau) <= A.b
    tband = np.abs(L - tau) <= A.b
    if tstar.sum() < 10 or tband.sum() < 10:
        return None
    inb = tstar & tband
    out = tstar & ~tband
    if inb.sum() < 5 or out.sum() < 5:
        return None
    T = max(int(A.frac * N), 1)
    lam = (A.sig ** 2) / (A.eta ** 2)
    idx_all = np.arange(N)
    idx_band = np.where(tband)[0]
    rows = {}
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
            s[j] += cnt[hit] * F[j] + A.sig * np.sqrt(cnt[hit]) * rs.standard_normal(j.size)
        seen = c > 0
        Fh = L.copy()
        Fh[seen] = (lam * L[seen] + s[seen]) / (lam + c[seen])
        wrong = (Fh >= tau) != H
        rows[al] = (float(wrong[tstar].sum()), float(wrong[tband].sum()),
                    float(wrong[inb].sum()), float(wrong[out].sum()))
    size = (float(tstar.sum()), float(tband.sum()), float(inb.sum()), float(out.sum()))
    return rows, size


acc = {al: np.zeros(4) for al in ALPHAS}
tot = np.zeros(4)
nseed = 0
t0 = time.time()
for sd in range(A.seeds):
    rs = np.random.default_rng(7717 + 31 * sd)
    got = run_one(rs)
    if got is None:
        continue
    r, size = got
    tot += np.array(size)
    nseed += 1
    for al in ALPHAS:
        acc[al] += np.array(r[al])

# Error rates are pooled over seeds, weighted by how many wells each seed
# contributes to each set, so the identity below holds exactly rather than up
# to the variation of the set sizes across seeds.
base = list(acc[0.0] / tot)
omega = float(tot[2] / tot[0])
res = {'setting': {'P': A.P, 'm': A.m, 'rank': A.rank, 'eta': A.eta, 'sig': A.sig,
                   'q': A.q, 'b': A.b, 'frac': A.frac, 'seeds': nseed},
       'omega': round(omega, 5), 'arms': {}}
for al in ALPHAS:
    e = list(acc[al] / tot)
    gS, gB, gI, gO = (base[0] - e[0], base[1] - e[1], base[2] - e[2], base[3] - e[3])
    kap = gI / gB if abs(gB) > 1e-9 else None
    mu = gO / gB if abs(gB) > 1e-9 else None
    pred = None
    if kap is not None and mu is not None:
        den = omega * kap + (1.0 - omega) * mu
        pred = 1.0 / den if abs(den) > 1e-9 else None
    res['arms'][str(al)] = {
        'errBstar': round(e[0], 5), 'errB': round(e[1], 5),
        'errIn': round(e[2], 5), 'errOut': round(e[3], 5),
        'gain_Bstar': round(gS, 5), 'gain_B': round(gB, 5),
        'gain_in': round(gI, 5), 'gain_out': round(gO, 5),
        'kappa': None if kap is None else round(kap, 4),
        'mu': None if mu is None else round(mu, 4),
        'factor_measured': round(gB / gS, 4) if abs(gS) > 2e-3 and gB > 0 else None,
        'factor_predicted': None if pred is None else round(pred, 4)}
    print(json.dumps({'k': 'tick', 'alpha': al, 'kappa': res['arms'][str(al)]['kappa'],
                      'mu': res['arms'][str(al)]['mu'],
                      'x': res['arms'][str(al)]['factor_measured'],
                      'xhat': res['arms'][str(al)]['factor_predicted']}), flush=True)

json.dump(res, open(A.out + '.tmp', 'w'), indent=1)
os.replace(A.out + '.tmp', A.out)
print(json.dumps({'k': 'done', 'out': A.out, 'sec': round(time.time() - t0, 1)}), flush=True)
