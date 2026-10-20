"""Intervals resampled over drug pairs, and the measured decomposition.

Two things are computed from the records written by wy.py.

First, intervals for the gain of a design over the non-targeting control on the
model-free band B*, the same gain on the model-defined band B, and the
difference between them, with the resampling unit being the drug pair rather
than the noise seed.  A screen can add or drop drug pairs; it cannot add noise
seeds, so an interval over pairs is the one that speaks to whether the effect
would survive a different panel.  Errors are summed within pairs by wy.py, so a
resample here means drawing pairs with replacement and recomputing both error
rates from the pooled counts.

Second, the three quantities of the ratio proposition, measured rather than
assumed: the overlap omega = |B n B*| / |B*|, and the two corrections kappa and
mu, defined as the gain on B n B* and the gain on B* \\ B, each over the gain on
B.  The predicted factor 1 / (omega kappa + (1 - omega) mu) is printed beside
the factor measured directly, which is a check on the algebra and on the code.

Usage:
    python3 pair_ci.py ww_oneil_Y_svd_homo.json --frac 0.32 --q 0.90 --b 0.30
"""
import argparse, json
import numpy as np

A = argparse.ArgumentParser()
A.add_argument('files', nargs='+')
A.add_argument('--frac', type=float, default=0.32)
A.add_argument('--q', type=float, default=0.90)
A.add_argument('--b', type=float, default=0.30)
A.add_argument('--ctrl', default='blind')
A.add_argument('--boot', type=int, default=4000)
A.add_argument('--seed', type=int, default=20260101)
A.add_argument('--arms', default='passive,band_hat,lse,ig,band_inf,ucb,straddle')
A.add_argument('--out', default='')
A = A.parse_args()

ARMS = A.arms.split(',')
LEVEL = (2.5, 97.5)


def close(x, y):
    return abs(float(x) - y) < 1e-9


def load(path):
    rows = json.load(open(path))
    den = [r for r in rows if r.get('meth') == 'pairdenom'
           and close(r['q'], A.q) and close(r['b'], A.b)]
    if not den:
        raise SystemExit('%s: no per-pair denominators at q=%g b=%g' % (path, A.q, A.b))
    sel = [r for r in rows if r.get('meth') != 'pairdenom'
           and close(r['q'], A.q) and close(r['b'], A.b) and close(r['frac'], A.frac)]
    if not sel:
        raise SystemExit('%s: no records at frac=%g' % (path, A.frac))
    return den[0], sel


def gather(sel, meth):
    """Per-pair error counts summed over seeds, for B* and for B."""
    rs = [r for r in sel if r['meth'] == meth and 'pe' in r]
    if not rs:
        return None
    pe = np.sum([np.asarray(r['pe'], float) for r in rs], axis=0)
    qe = np.sum([np.asarray(r['qe'], float) for r in rs], axis=0)
    return pe, qe, len(rs)


def pct(v):
    return [round(float(np.percentile(v, LEVEL[0])), 5),
            round(float(np.percentile(v, LEVEL[1])), 5)]


out = []
for path in A.files:
    den, sel = load(path)
    nT = np.asarray(den['pnT'], float)
    nB = np.asarray(den['pnB'], float)
    P = nT.size
    ctrl = gather(sel, A.ctrl)
    if ctrl is None:
        raise SystemExit('%s: control arm %s carries no per-pair counts' % (path, A.ctrl))
    cpe, cqe, nseed = ctrl
    rng = np.random.default_rng(A.seed)
    draw = rng.integers(0, P, size=(A.boot, P))
    # The resampled denominators are shared by every arm, so build them once.
    dT = nT[draw].sum(axis=1) * nseed
    dB = nB[draw].sum(axis=1) * nseed
    cET = cpe[draw].sum(axis=1)
    cEB = cqe[draw].sum(axis=1)
    ds = sel[0]['ds']
    base = {'file': path, 'ds': ds, 'q': A.q, 'b': A.b, 'frac': A.frac,
            'pairs': int(P), 'seeds': nseed, 'boot': A.boot,
            'MT': sel[0]['MT'], 'M': sel[0]['M'], 'MI': sel[0]['MI'], 'MO': sel[0]['MO'],
            'omega': round(float(sel[0]['MI']) / float(sel[0]['MT']), 5)}
    omega = base['omega']
    rows = []
    for mt in ARMS:
        got = gather(sel, mt)
        if got is None:
            continue
        ape, aqe, _ = got
        gT = (cET - ape[draw].sum(axis=1)) / dT
        gB = (cEB - aqe[draw].sum(axis=1)) / dB
        pointT = float((cpe.sum() - ape.sum()) / (nT.sum() * nseed))
        pointB = float((cqe.sum() - aqe.sum()) / (nB.sum() * nseed))
        diff = gB - gT
        ok = gT > 0
        ratio = np.where(ok, gB / np.where(ok, gT, 1.0), np.nan)
        frac_ok = float(np.mean(ok))
        ectrl = [r for r in sel if r['meth'] == A.ctrl]
        earm = [r for r in sel if r['meth'] == mt]
        gi = float(np.mean([r['errI'] for r in ectrl]) - np.mean([r['errI'] for r in earm]))
        go = float(np.mean([r['errO'] for r in ectrl]) - np.mean([r['errO'] for r in earm]))
        gb = float(np.mean([r['berr'] for r in ectrl]) - np.mean([r['berr'] for r in earm]))
        kap = gi / gb if abs(gb) > 1e-9 else None
        mu = go / gb if abs(gb) > 1e-9 else None
        pred = None
        if kap is not None and mu is not None:
            d = omega * kap + (1.0 - omega) * mu
            pred = round(1.0 / d, 4) if abs(d) > 1e-9 else None
        row = dict(base)
        row.update({
            'meth': mt,
            'gain_Bstar': round(pointT, 5), 'gain_Bstar_ci': pct(gT),
            'gain_B': round(pointB, 5), 'gain_B_ci': pct(gB),
            'diff': round(pointB - pointT, 5), 'diff_ci': pct(diff),
            'diff_excludes_zero': bool(np.percentile(diff, LEVEL[0]) > 0),
            'gain_in': round(gi, 5), 'gain_out': round(go, 5), 'gain_band': round(gb, 5),
            'kappa': None if kap is None else round(kap, 4),
            'mu': None if mu is None else round(mu, 4),
            'factor_measured': round(pointB / pointT, 4) if pointT > 2e-4 else None,
            'factor_predicted': pred,
            'factor_ci': pct(ratio[ok]) if frac_ok >= 0.90 else None,
            'frac_usable': round(frac_ok, 4)})
        rows.append(row)
        print(json.dumps({k: row[k] for k in
                          ('ds', 'meth', 'gain_Bstar', 'gain_Bstar_ci', 'gain_B',
                           'gain_B_ci', 'diff', 'diff_ci', 'omega', 'kappa', 'mu',
                           'factor_measured', 'factor_predicted', 'factor_ci')}),
              flush=True)
    out.extend(rows)

if A.out:
    json.dump(out, open(A.out, 'w'), indent=1)
    print(json.dumps({'k': 'done', 'rows': len(out), 'out': A.out}), flush=True)
