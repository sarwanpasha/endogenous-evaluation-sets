"""Build zy.py from xy.py.

Three additions, all requested by the review:

  1. berrH -- boundary error on a band built from the stage-one fit Lh, which
     uses only purchased data, instead of the reference completion Lv, which is
     a truncated SVD of the complete matrix and therefore sees the values the
     designs are scored against.
  2. hit-rate at k -- recall of the true top-k wells by the fused score, which
     is the endpoint the deployed adaptive screens actually report.
  3. cost-aware arms -- the Bayes rule for an asymmetric cost k thresholds the
     posterior probability at 1/(1+k) rather than the fused mean at tau.

Everything else, including the random seed stream, is untouched, so the
existing arms reproduce the published numbers exactly.
"""
import re, sys, os

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'xy.py')
DST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zy.py')
s = open(SRC).read()
fails = []


def sub(pat, rep, n=1, flags=0):
    global s
    s2, k = re.subn(pat, rep, s, count=n, flags=flags)
    if k != n:
        fails.append(pat[:70])
    s = s2


# ---------------------------------------------------------------- metrics ---
sub(r'def met\(pred, Ht, bm, tm, bal\):',
    'def met(pred, Ht, bm, tm, bal, hm=None, sc=None):')

sub(r"(\n)(    )return o\n",
    r"\1\2if hm is not None and bool(hm.any()):\n"
    r"\2    o['berrH'] = round(float((pred[hm] != Ht[hm]).mean()), 5)\n"
    r"\2    o['MH'] = int(hm.sum())\n"
    r"\2if sc is not None:\n"
    r"\2    npos = int(Ht.sum())\n"
    r"\2    order = np.argsort(-sc, kind='stable')\n"
    r"\2    for kk in (100, 500, npos):\n"
    r"\2        if kk <= 0 or kk > sc.size: continue\n"
    r"\2        nm = 'hitpos' if kk == npos else 'hit%d' % kk\n"
    r"\2        o[nm] = round(float(Ht[order[:kk]].mean()), 5)\n"
    r"\1\2return o\n")

# ------------------------------------------------- posterior tail probability -
sub(r"(\n# ---------------------------------------------------------------- metrics\n)",
    "\n\ndef postp(prior, s, cnt, sig, sdv, tau):\n"
    "    \"\"\"Posterior probability that the entry is above tau, under the same\n"
    "    normal-normal model the point rule uses.\"\"\"\n"
    "    ok = cnt > 0\n"
    "    yb = np.zeros(cnt.size); yb[ok] = s[ok] / cnt[ok]\n"
    "    if int(ok.sum()) > 20:\n"
    "        rr = float(np.mean((yb[ok] - prior[ok]) ** 2)\n"
    "                   - sig ** 2 * float(np.mean(1.0 / cnt[ok])))\n"
    "        vh = max(rr, 1e-4 * sdv ** 2)\n"
    "    else:\n"
    "        vh = sdv ** 2\n"
    "    prec = 1.0 / vh + cnt / sig ** 2\n"
    "    mp = (prior / vh + s / sig ** 2) / prec\n"
    "    return _nrm.cdf((mp - tau) * np.sqrt(prec))\n"
    r"\1")

# ------------------------------------------------------- the evaluation band --
sub(r"(\n                )tm = np\.abs\(Fv - tau\) <= b \* sd\n",
    r"\1tm = np.abs(Fv - tau) <= b * sd"
    r"\1hm = np.abs(Lh - tau) <= b * sd\n")

# -------------------------------------------------------------- recording ----
sub(r"def rec\(mt, pred, extra=None\):", 'def rec(mt, pred, extra=None, sc=None):')
sub(r"row\.update\(met\(pred, Ht, bm, tm, bal\)\)",
    'row.update(met(pred, Ht, bm, tm, bal, hm, sc))')

# ------------------------------------------------- capture the fused scores ---
sub(r"vb = post\(MEANP, sm_, cn, SIGSC, sd\) >= tau",
    'svb = post(MEANP, sm_, cn, SIGSC, sd); vb = svb >= tau')
sub(r"vp = post\(Lp, sm_, cn, SIGSC, sd\) >= tau",
    'svp = post(Lp, sm_, cn, SIGSC, sd); vp = svp >= tau')
sub(r"rec\('blind', vb, \{", "rec('blind', vb, sc=svb, extra={")
sub(r"rec\('passive', vp, \{", "rec('passive', vp, sc=svp, extra={")
sub(r"rec\('band_inf', post\(Lv, si, ci, SIGSC, sd\) >= tau\)",
    "sbi = post(Lv, si, ci, SIGSC, sd)\n"
    "                rec('band_inf', sbi >= tau, sc=sbi)")
sub(r"rec\('band_hat', post\(Lh, st, ct, SIGSC, sd\) >= tau, \{'Mhat': int\(arh\.size\)\}\)",
    "sbh = post(Lh, st, ct, SIGSC, sd)\n"
    "                rec('band_hat', sbh >= tau, {'Mhat': int(arh.size)}, sc=sbh)")
sub(r"rec\('lse', post\(Lh, st, ct, SIGSC, sd\) >= tau, \{'Mhat': int\(arh\.size\)\}\)",
    "sls = post(Lh, st, ct, SIGSC, sd)\n"
    "                rec('lse', sls >= tau, {'Mhat': int(arh.size)}, sc=sls)")

# ------------------------------------------------------- cost-aware arms ------
sub(r"(\n                if DESIGNS is None or 'ig' in DESIGNS:)",
    "\n                for _k in (3, 10):\n"
    "                    _thr = 1.0 / (1.0 + _k)\n"
    "                    rec('blind_k%d' % _k,\n"
    "                        postp(MEANP, sm_, cn, SIGSC, sd, tau) >= _thr)\n"
    "                    rec('passive_k%d' % _k,\n"
    "                        postp(Lp, sm_, cn, SIGSC, sd, tau) >= _thr)\n"
    "                    rec('band_inf_k%d' % _k,\n"
    "                        postp(Lv, si, ci, SIGSC, sd, tau) >= _thr)\n"
    "                    rec('band_hat_k%d' % _k,\n"
    "                        postp(Lh, st, ct, SIGSC, sd, tau) >= _thr)"
    r"\1")

# ------------------------------------------------------------- output name ----
sub(r"OUTF = 'xy_", "OUTF = 'zz_")

open(DST, 'w').write(s)
if fails:
    print('FAILED:', fails)
    sys.exit(1)
import ast
ast.parse(s)
print('wrote %s, %d lines, parses clean' % (DST, s.count('\n') + 1))
