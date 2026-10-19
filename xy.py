"""Extended driver: adds trivial baselines, decision-cost metrics, a class-balanced
boundary set, an empirical heteroskedastic noise model, a band-localised
classification prior, and an information-gain adaptive design."""
import numpy as np, sys, json, os, time, argparse
exec(open('x12.py').read().split(chr(10) + "DS=[('oneil")[0])
import feat
from scipy.stats import norm as _nrm

def pred_err(prior, cnt, tau, mask, sig, lam):
    D = Fv[mask] - tau; pm = prior[mask]; dl = pm - Fv[mask]
    c = cnt[mask].astype(float); sg = np.where(D >= 0, 1.0, -1.0)
    out = np.empty(D.size); z = c <= 0
    if z.any():
        out[z] = ((pm[z] >= tau) != (Fv[mask][z] >= tau)).astype(float)
    nz = ~z
    if nz.any():
        w = lam / (lam + c[nz])
        arg = -(lam + c[nz]) / (sig * np.sqrt(c[nz])) * (np.abs(D[nz]) + w * dl[nz] * sg[nz])
        out[nz] = _nrm.cdf(arg)
    return float(out.mean())

def lam_of(prior, s, cnt, sig, sd):
    ok = cnt > 0
    if int(ok.sum()) > 20:
        yb = s[ok] / cnt[ok]
        rr = float(np.mean((yb - prior[ok]) ** 2) - sig ** 2 * float(np.mean(1.0 / cnt[ok])))
        vh = max(rr, 1e-4 * sd ** 2)
    else:
        vh = sd ** 2
    return sig ** 2 / vh

ap = argparse.ArgumentParser()
ap.add_argument('--ds', required=True)
ap.add_argument('--score', default='')
ap.add_argument('--prior', default='svd')          # svd | gbt | mlp | clf | clfu
ap.add_argument('--seeds', type=int, default=10)
ap.add_argument('--sigf', type=float, default=0.50)
ap.add_argument('--noise', default='homo')          # homo | emp
ap.add_argument('--rank', type=int, default=0)
ap.add_argument('--fracs', default='0.01,0.04,0.16,0.32,1.28,5.12')
ap.add_argument('--bands', default='0.10,0.30')
ap.add_argument('--quants', default='0.90')
ap.add_argument('--a1', type=float, default=0.30)
ap.add_argument('--kappa', type=float, default=0.5)   # band width of the clf weight
ap.add_argument('--gamma', type=float, default=1.0)   # clf prior offset, in sd units
ap.add_argument('--igrounds', type=int, default=120)
ap.add_argument('--designs', default='all')
ap.add_argument('--tag', default='')
A = ap.parse_args()

CFG = {'oneil': ('df_OncoPolyPharmacology.pkl', 'Y', 8), 'dcomb': ('df_DrugComb.pkl', 'Synergy_Loewe', 12)}
pth, sdef, rdef = CFG[A.ds]
scol = A.score if A.score else sdef
r = A.rank if A.rank > 0 else rdef
FR = [float(x) for x in A.fracs.split(',')]
BW = [float(x) for x in A.bands.split(',')]
QS = [float(x) for x in A.quants.split(',')]
DESIGNS = None if A.designs == 'all' else set(A.designs.split(','))

t0 = time.time()
F, keys, cells, SM, CE = feat.load(pth, scol)
vm = np.isfinite(F); Fv = F[vm]; N = int(vm.sum()); P, m = F.shape
sd = float(Fv.std()); sbar = (A.sigf if A.sigf > 0 else 0.50) * sd
ix = np.argwhere(vm); i0 = ix[:, 0]; j0 = ix[:, 1]; K = N
G = F.copy(); G[~vm] = Fv.mean()
Uu, Ss, Vv = np.linalg.svd(G, full_matrices=False)
Lv = ((Uu[:, :r] * Ss[:r]) @ Vv[:r])[vm]
TAU = {q: float(np.quantile(Fv, q)) for q in QS}
BAND = {(q, b): np.where(np.abs(Lv - TAU[q]) <= b * sd) [0] for q in QS for b in BW}

# ---------------------------------------------------------------- noise model
def emp_sigma():
    """Per-entry sd from the screen's own repeated wells.  We bin the repeated
    wells by |mean| and use the observed within-well sd in each bin, which is
    model free and so does not impose monotonicity; entries are then assigned by
    |F_e|.  If --sigf is 0 the absolute calibrated level is used, otherwise the
    curve is rescaled so that the mean variance matches sigf*sd."""
    import pandas as pd
    d = pd.read_pickle(pth)
    a = d['Drug1_ID'].astype(str).values; bb = d['Drug2_ID'].astype(str).values
    kk = np.where(a < bb, a + '|' + bb, bb + '|' + a)
    t = pd.DataFrame({'k': kk, 'c': d['Cell_Line_ID'].astype(str).values,
                      'y': pd.to_numeric(d[scol], errors='coerce').values}).dropna()
    gp = t.groupby(['k', 'c'])['y']; n = gp.size(); mu = gp.mean(); sv = gp.std(ddof=1)
    sel = (n >= 2) & sv.notna()
    info = {'groups_with_reps': int(sel.sum()), 'frac_with_reps': round(float((n >= 2).mean()), 4)}
    del d, t
    if int(sel.sum()) < 50:
        info['fallback'] = 'homoskedastic (too few repeated wells)'
        return np.full(N, sbar), info
    nn = n[sel].astype(float).values; sw = sv[sel].values; mm = np.abs(mu[sel].values)
    pooled = float(np.sqrt(np.sum((nn - 1) * sw ** 2) / np.sum(nn - 1)))
    info['pooled_within_sd'] = round(pooled, 4)
    info['calibrated_sigma_over_sd'] = round(pooled / sd, 4)
    nb = 10
    edges = np.quantile(mm, np.linspace(0, 1, nb + 1))
    bidx = np.clip(np.digitize(mm, edges[1:-1]), 0, nb - 1)
    curve = np.array([np.sqrt(np.mean(sw[bidx == j] ** 2)) if (bidx == j).sum() > 5 else pooled
                      for j in range(nb)])
    info['sd_by_absmean_decile'] = [round(float(z), 4) for z in curve]
    info['decile_ratio'] = round(float(curve.max() / max(curve.min(), 1e-9)), 3)
    g = curve[np.clip(np.digitize(np.abs(Fv), edges[1:-1]), 0, nb - 1)]
    target = pooled if A.sigf <= 0 else sbar
    s_e = g * target / np.sqrt(np.mean(g ** 2))
    info['target_rms_sigma'] = round(float(target), 4)
    info['target_over_sd'] = round(float(target / sd), 4)
    info['sigma_p90_over_p10'] = round(float(np.quantile(s_e, .9) / np.quantile(s_e, .1)), 3)
    return s_e, info

NOISE_INFO = {}
if A.noise == 'emp':
    SIG, NOISE_INFO = emp_sigma()
else:
    SIG = np.full(N, sbar)
SIGSC = float(np.sqrt(np.mean(SIG ** 2)))     # scalar used by the estimator

def draw2(rs, cnt, fv, sg):
    s = np.zeros(cnt.size); ok = cnt > 0
    s[ok] = cnt[ok] * fv[ok] + sg[ok] * np.sqrt(cnt[ok]) * rs.normal(0, 1, int(ok.sum()))
    return s

# ---------------------------------------------------------------- features
def build_feats():
    dset = sorted({d for k in keys for d in str(k).split('|')})
    fp = {}; mode = 'onehot'
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFingerprintGenerator as FG
        gen = FG.GetMorganGenerator(radius=2, fpSize=1024); ok = 0
        for dd in dset:
            s = SM.get(dd, ''); mo = Chem.MolFromSmiles(s) if s else None
            if mo is None: fp[dd] = np.zeros(1024, np.float32)
            else: fp[dd] = gen.GetFingerprintAsNumPy(mo).astype(np.float32); ok += 1
        if ok < 0.5 * len(dset): raise RuntimeError('too few parsed')
        mode = 'morgan1024_r2_parsed%d/%d' % (ok, len(dset))
    except Exception as e:
        di = {d: i for i, d in enumerate(dset)}
        for dd in dset:
            z = np.zeros(len(dset), np.float32); z[di[dd]] = 1.0; fp[dd] = z
        mode = 'onehot%d (%s)' % (len(dset), str(e)[:40])
    Xp = np.stack([fp[str(k).split('|')[0]] + fp[str(k).split('|')[1]] for k in keys]).astype(np.float32)
    if CE:
        E = np.stack([CE[str(c)] for c in cells]).astype(np.float32)
        from sklearn.decomposition import PCA
        nq = int(min(32, E.shape[0] - 1, E.shape[1]))
        Xc = PCA(n_components=nq, random_state=0).fit_transform(E - E.mean(0)).astype(np.float32)
        cm = 'expPCA%d_from%d' % (nq, E.shape[1])
    else:
        Xc = np.eye(len(cells), dtype=np.float32); cm = 'cellonehot%d' % len(cells)
    Xc = (Xc - Xc.mean(0)) / (Xc.std(0) + 1e-8)
    return np.hstack([Xp[i0], Xc[j0]]).astype(np.float32), mode, cm

X = None; FMODE = CMODE = 'n/a'
if A.prior in ('gbt', 'mlp', 'clf', 'clfu'):
    X, FMODE, CMODE = build_feats()

# ---------------------------------------------------------------- priors
def prior_spine(cn, sm_):
    Ym = np.zeros((P, m)); Cc = np.zeros((P, m))
    Ym[vm] = np.where(cn > 0, sm_ / np.maximum(cn, 1), 0.0); Cc[vm] = cn
    return fit(Ym, Cc, r)[vm]

def prior_reg(cn, sm_, seed):
    q = cn > 0; yb = (sm_[q] / cn[q]).astype(np.float32)
    from sklearn.ensemble import HistGradientBoostingRegressor
    md = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, random_state=seed)
    md.fit(X[q], yb, sample_weight=cn[q].astype(np.float64))
    return md.predict(X).astype(float)

def prior_clf(cn, sm_, seed, tau, localise):
    """Train for the SIGN of F-tau near tau, not for squared error.  The prior is
    placed at tau +/- gamma*sd according to the predicted class probability, so
    sign(L-tau) is exactly the classifier's decision."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    q = cn > 0
    yb = sm_[q] / cn[q]
    z = (yb >= tau).astype(int)
    w = cn[q].astype(np.float64)
    if localise:
        w = w * np.exp(-0.5 * ((yb - tau) / (A.kappa * sd)) ** 2)
    if z.min() == z.max() or w.sum() <= 0:
        return np.full(N, float(Fv.mean()))
    md = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.08, random_state=seed)
    md.fit(X[q], z, sample_weight=w)
    p = md.predict_proba(X)[:, 1].astype(float)
    return tau + A.gamma * sd * (2.0 * p - 1.0)

def prior_fit(cn, sm_, seed, tau=None):
    if A.prior == 'svd': return prior_spine(cn, sm_)
    if A.prior == 'gbt': return prior_reg(cn, sm_, seed)
    if A.prior in ('clf', 'clfu'): return prior_clf(cn, sm_, seed, tau, A.prior == 'clf')
    raise SystemExit('prior %s not supported here' % A.prior)

# ---------------------------------------------------------------- info gain
_GH_X, _GH_W = np.polynomial.hermite_e.hermegauss(7)
_GH_W = _GH_W / _GH_W.sum()

def _H(p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(p * np.log(p) + (1 - p) * np.log1p(-p))

def ig_design(rs, tau, Lh, vh, T, c0, s0, rounds):
    """Batched design maximising expected reduction in the entropy of the
    thresholded label 1{F>=tau} under the working posterior."""
    cnt = c0.copy().astype(np.int64); s = s0.copy()
    rem = int(T - cnt.sum())
    if rem <= 0: return s, cnt
    bs = max(1, rem // max(rounds, 1))
    v = max(vh, 1e-8 * sd ** 2)
    while rem > 0:
        b = int(min(bs, rem))
        prec = 1.0 / v + cnt / SIGSC ** 2
        mpost = (Lh / v + s / SIGSC ** 2) / prec
        t = 1.0 / np.sqrt(prec)
        p = _nrm.cdf((mpost - tau) / t)
        h0 = _H(p)
        sy = np.sqrt(t ** 2 + SIGSC ** 2)          # predictive sd of one new draw
        pr2 = prec + 1.0 / SIGSC ** 2
        t2 = 1.0 / np.sqrt(pr2)
        eh = np.zeros(N)
        for gx, gw in zip(_GH_X, _GH_W):
            y = mpost + sy * gx
            m2 = (mpost * prec + y / SIGSC ** 2) / pr2
            eh += gw * _H(_nrm.cdf((m2 - tau) / t2))
        gain = h0 - eh
        sel = np.argpartition(-gain, b - 1)[:b] if b < N else np.arange(N)
        s[sel] += Fv[sel] + SIG[sel] * rs.normal(0, 1, sel.size)
        cnt[sel] += 1; rem -= b
    return s, cnt

# ---------------------------------------------------------------- metrics
def met(pred, Ht, bm, tm, bal):
    tp = float((pred & Ht).sum()); fp = float((pred & ~Ht).sum()); fn = float(((~pred) & Ht).sum())
    f1v = 0.0 if tp == 0 else 2 * tp / (2 * tp + fp + fn)
    pr = 0.0 if tp + fp == 0 else tp / (tp + fp)
    rc = 0.0 if tp + fn == 0 else tp / (tp + fn)
    pt = pred[tm]; ht = Ht[tm]
    btp = float((pt & ht).sum()); bfp = float((pt & ~ht).sum()); bfn = float(((~pt) & ht).sum())
    nb = float(tm.sum())
    o = {'f1': round(f1v, 5), 'prec': round(pr, 5), 'rec': round(rc, 5),
         'gerr': round(float((pred != Ht).mean()), 5),
         'gfp': round(fp / pred.size, 6), 'gfn': round(fn / pred.size, 6),
         'nsel': round(float(pred.mean()), 6),
         'berr': round(float((pred[bm] != Ht[bm]).mean()), 5),
         'berrT': round(float((pt != ht).mean()), 5),
         'bfp': round(bfp / nb, 5), 'bfn': round(bfn / nb, 5),
         'bprec': round(0.0 if btp + bfp == 0 else btp / (btp + bfp), 5),
         'brec': round(0.0 if btp + bfn == 0 else btp / (btp + bfn), 5)}
    for kk in (3, 10):
        o['cost%d' % kk] = round((bfp + kk * bfn) / nb, 5)
        o['gcost%d' % kk] = round((fp + kk * fn) / pred.size, 6)
    if bal is not None and bal.size:
        o['berrBal'] = round(float((pred[bal] != Ht[bal]).mean()), 5)
    return o

print(json.dumps({'k': 'meta', 'ds': A.ds, 'score': scol, 'prior': A.prior, 'noise': A.noise,
                  'N': N, 'P': P, 'm': m, 'r': r, 'sd': sd, 'sigf': A.sigf, 'sbar': sbar,
                  'sigsc': SIGSC, 'noise_info': NOISE_INFO, 'seeds': A.seeds, 'a1': A.a1,
                  'kappa': A.kappa, 'gamma': A.gamma, 'feat': FMODE, 'cellfeat': CMODE,
                  'fracs': FR, 'bands': BW, 'quants': QS,
                  'load_sec': round(time.time() - t0, 1)}), flush=True)

OUTF = 'xy_%s_%s_%s_%s%s.json' % (A.ds, scol, A.prior, A.noise, A.tag)
res = []; done = set()
if os.path.exists(OUTF):
    try:
        res = json.load(open(OUTF))
        done = {(round(float(x['frac']), 6), int(x['seed'])) for x in res}
        print(json.dumps({'k': 'resume', 'rows': len(res), 'cells': len(done)}), flush=True)
    except Exception as e:
        print(json.dumps({'k': 'resume_failed', 'err': str(e)[:80]}), flush=True); res = []; done = set()

# class-balanced boundary subsets, fixed across seeds
BAL = {}
for q in QS:
    for b in BW:
        tmk = np.where(np.abs(Fv - TAU[q]) <= b * sd)[0]
        hp = tmk[Fv[tmk] >= TAU[q]]; hn = tmk[Fv[tmk] < TAU[q]]
        n2 = int(min(hp.size, hn.size))
        if n2 > 0:
            g = np.random.default_rng(12345)
            BAL[(q, b)] = np.concatenate([g.choice(hp, n2, replace=False), g.choice(hn, n2, replace=False)])
        else:
            BAL[(q, b)] = np.array([], dtype=int)

for frac in FR:
    T = max(int(frac * N), 1)
    for sdd in range(A.seeds):
        if (round(float(frac), 6), int(sdd)) in done: continue
        ts = time.time()
        rs = np.random.default_rng(8317 * (int(frac * 100000) + 1) + sdd * 7919 + int(A.sigf * 1000))
        cn = np.bincount(rs.integers(0, K, T), minlength=K); sm_ = draw2(rs, cn, Fv, SIG)
        T1 = max(int(A.a1 * T), 1)
        c1 = np.bincount(rs.integers(0, K, T1), minlength=K); s1 = draw2(rs, c1, Fv, SIG)
        MEANP = np.full(K, float(Fv.mean()))
        for q in QS:
            tau = TAU[q]; Ht = Fv >= tau
            Lp = prior_fit(cn, sm_, sdd, tau)
            Lh = prior_fit(c1, s1, 1000 + sdd, tau)
            lamB = lam_of(MEANP, sm_, cn, SIGSC, sd)
            lamP = lam_of(Lp, sm_, cn, SIGSC, sd)
            vhat = SIGSC ** 2 / max(lam_of(Lh, s1, c1, SIGSC, sd), 1e-9)
            vb = post(MEANP, sm_, cn, SIGSC, sd) >= tau
            vp = post(Lp, sm_, cn, SIGSC, sd) >= tau
            vmeas = np.zeros(K, dtype=bool); oq = cn > 0
            vmeas[oq] = (sm_[oq] / cn[oq]) >= tau
            for b in BW:
                ar = BAND[(q, b)]
                if ar.size == 0: continue
                bm = np.zeros(K, dtype=bool); bm[ar] = True
                tm = np.abs(Fv - tau) <= b * sd
                MT = int(tm.sum())
                if MT == 0: continue
                bal = BAL[(q, b)]
                base = {'ds': A.ds, 'score': scol, 'prior': A.prior, 'noise': A.noise,
                        'q': q, 'b': b, 'frac': frac, 'seed': sdd, 'M': int(ar.size), 'MT': MT,
                        'posfrac': round(float(Ht[tm].mean()), 5),
                        'pri_e': round(float(((Lv >= tau) != Ht)[tm].mean()), 5),
                        'prip_e': round(float(((Lp >= tau) != Ht)[tm].mean()), 5),
                        'prih_e': round(float(((Lh >= tau) != Ht)[tm].mean()), 5)}
                def rec(mt, pred, extra=None):
                    if DESIGNS is not None and mt not in DESIGNS: return
                    row = dict(base); row['meth'] = mt
                    row.update(met(pred, Ht, bm, tm, bal))
                    if extra: row.update(extra)
                    res.append(row)
                rec('null0', np.zeros(K, dtype=bool))
                rec('null1', np.ones(K, dtype=bool))
                rec('meas', vmeas)
                rec('blind', vb, {'lam': round(lamB, 4), 'pred': round(pred_err(MEANP, cn, tau, tm, SIGSC, lamB), 5)})
                rec('passive', vp, {'lam': round(lamP, 4), 'pred': round(pred_err(Lp, cn, tau, tm, SIGSC, lamP), 5)})
                ci = np.zeros(K, dtype=np.int64); si = np.zeros(K); c2 = alloc(rs, ar.size, T)
                ci[ar] = c2; si[ar] = draw2(rs, c2, Fv[ar], SIG[ar])
                rec('band_inf', post(Lv, si, ci, SIGSC, sd) >= tau)
                arh = np.where(np.abs(Lh - tau) <= 2 * b * sd)[0]
                rem = T - T1
                ct = c1.copy(); st = s1.copy()
                if arh.size > 0 and rem > 0:
                    c3 = alloc(rs, arh.size, rem); ct[arh] += c3; st[arh] += draw2(rs, c3, Fv[arh], SIG[arh])
                rec('band_hat', post(Lh, st, ct, SIGSC, sd) >= tau, {'Mhat': int(arh.size)})
                ct = c1.copy(); st = s1.copy()
                if arh.size > 0 and rem > 0:
                    s2, c2b = apt2(rs, Fv[arh], SIGSC, rem, tau, 0.02 * sd); ct[arh] += c2b; st[arh] += s2
                rec('lse', post(Lh, st, ct, SIGSC, sd) >= tau, {'Mhat': int(arh.size)})
                # sample-split arm: prior from stage one, fused ONLY with the
                # disjoint stage-two draws, so L is independent of the
                # measurements it is fused with and Thm 1 applies exactly.
                cu = np.zeros(K, dtype=np.int64); su = np.zeros(K)
                if rem > 0:
                    cu = np.bincount(rs.integers(0, K, rem), minlength=K)
                    su = draw2(rs, cu, Fv, SIG)
                lamS = lam_of(Lh, su, cu, SIGSC, sd)
                rec('passive_split', post(Lh, su, cu, SIGSC, sd) >= tau,
                    {'lam': round(lamS, 4),
                     'pred': round(pred_err(Lh, cu, tau, tm, SIGSC, lamS), 5)})
                rec('blind_split', post(MEANP, su, cu, SIGSC, sd) >= tau,
                    {'pred': round(pred_err(MEANP, cu, tau, tm, SIGSC,
                                            lam_of(MEANP, su, cu, SIGSC, sd)), 5)})
                if DESIGNS is None or 'ig' in DESIGNS:
                    sg, cg = ig_design(rs, tau, Lh, vhat, T, c1, s1, A.igrounds)
                    rec('ig', post(Lh, sg, cg, SIGSC, sd) >= tau)
        print(json.dumps({'k': 'tick', 'frac': frac, 'seed': sdd, 'sec': round(time.time() - ts, 1)}), flush=True)
    json.dump(res, open(OUTF + '.tmp', 'w')); os.replace(OUTF + '.tmp', OUTF)
print(json.dumps({'k': 'done', 'rows': len(res), 'out': OUTF, 'sec': round(time.time() - t0, 1)}), flush=True)
