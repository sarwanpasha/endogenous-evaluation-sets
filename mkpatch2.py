"""Build wy.py from zy.py.

Four additions, all answering specific points raised in review:

  1. Subset error rates.  Every arm is scored separately on the overlap
     B n B*, on the part of the model-free band the model-defined band misses,
     B* \\ B, and on B itself.  With the non-targeting control these give the
     three quantities omega, kappa and mu that the ratio proposition names, so
     the decomposition is measured rather than assumed.
  2. Per-pair error counts.  Errors on B* and on B are summed within each drug
     pair, which lets an interval be resampled over drug pairs instead of over
     noise seeds.  Pairs, not wells, are the unit a screen can add or drop.
  3. A upper-confidence-bound design, which spends the second stage where the
     posterior mean plus two standard deviations is largest.  This is the
     exploitation loop that deployed recommenders run.
  4. A straddle design, which spends the second stage where 1.96 s(e) - |m(e) -
     tau| is largest.  This is the level-set acquisition of the active
     level-set estimation literature, applied to the same posterior.

The two new designs draw from their own generator, so the arms that exist in
zy.py reproduce identical numbers when wy.py is run on the same lattice.  The
per-pair vectors are long, so they are written only at one band and one budget,
set by --pairband and --pairfrac, which leaves the lattice itself free to match
the lattice the earlier runs used.
"""
import re, sys, os

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zy.py')
DST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'wy.py')
s = open(SRC).read()
fails = []


def sub(pat, rep, n=1, flags=0):
    global s
    s2, k = re.subn(pat, rep, s, count=n, flags=flags)
    if k != n:
        fails.append(pat[:70])
    s = s2


# ------------------------------------------------------------------ options --
sub(r"ap\.add_argument\('--tag', default=''\)",
    "ap.add_argument('--tag', default='')\n"
    "ap.add_argument('--pairs', type=int, default=1)\n"
    "ap.add_argument('--pairband', type=float, default=0.30)\n"
    "ap.add_argument('--pairfrac', type=float, default=0.32)")

sub(r"(DESIGNS = None if A\.designs == 'all' else set\(A\.designs\.split\(','\)\)\n)",
    r"\1PAIRW = A.pairs > 0\n")

# ------------------------------------- per-pair bookkeeping and two designs ---
sub(r"(\ndef postp\(prior, s, cnt, sig, sdv, tau\):)",
    "\n\nPAIRM = ('blind', 'passive', 'band_hat', 'lse', 'ig', 'band_inf', 'ucb',\n"
    "         'straddle')\n"
    "\n"
    "def pv(pred, Ht, msk):\n"
    "    \"\"\"Errors on the masked wells, summed within each drug pair.\"\"\"\n"
    "    w = (pred[msk] != Ht[msk]).astype(float)\n"
    "    return np.bincount(i0[msk], weights=w, minlength=P).astype(np.int64).tolist()\n"
    "\n"
    "def pn(msk):\n"
    "    \"\"\"How many of the masked wells each drug pair contributes.\"\"\"\n"
    "    return np.bincount(i0[msk], minlength=P).astype(np.int64).tolist()\n"
    "\n"
    "def greedy(rs, tau, Lh, vh, T, c0, s0, rounds, kind):\n"
    "    \"\"\"Batched greedy second stage under the working posterior.  kind 'ucb'\n"
    "    scores by the posterior mean plus two standard deviations; kind\n"
    "    'straddle' scores by 1.96 s(e) - |m(e) - tau|.\"\"\"\n"
    "    cnt = c0.copy().astype(np.int64); s = s0.copy()\n"
    "    rem = int(T - cnt.sum())\n"
    "    if rem <= 0: return s, cnt\n"
    "    bs = max(1, rem // max(rounds, 1))\n"
    "    v = max(vh, 1e-8 * sd ** 2)\n"
    "    while rem > 0:\n"
    "        b = int(min(bs, rem))\n"
    "        prec = 1.0 / v + cnt / SIGSC ** 2\n"
    "        mpost = (Lh / v + s / SIGSC ** 2) / prec\n"
    "        t = 1.0 / np.sqrt(prec)\n"
    "        sc = mpost + 2.0 * t if kind == 'ucb' else 1.96 * t - np.abs(mpost - tau)\n"
    "        sel = np.argpartition(-sc, b - 1)[:b] if b < N else np.arange(N)\n"
    "        s[sel] += Fv[sel] + SIG[sel] * rs.normal(0, 1, sel.size)\n"
    "        cnt[sel] += 1; rem -= b\n"
    "    return s, cnt\n"
    r"\1")

# ------------------------------------------------------------ subset errors --
sub(r"(\n)(    )return o\n",
    r"\1\2ii = tm & bm\n"
    r"\2if bool(ii.any()):\n"
    r"\2    o['errI'] = round(float((pred[ii] != Ht[ii]).mean()), 5)\n"
    r"\2oo = tm & ~bm\n"
    r"\2if bool(oo.any()):\n"
    r"\2    o['errO'] = round(float((pred[oo] != Ht[oo]).mean()), 5)\n"
    r"\1\2return o\n")

# --------------------------------------------------------- the denominators --
sub(r"res = \[\]; done = set\(\)\n", "res = []; done = set(); PDONE = set()\n")

sub(r"(\n                )bal = BAL\[\(q, b\)\]\n",
    r"\1bal = BAL[(q, b)]"
    r"\1MI = int((tm & bm).sum()); MO = int((tm & ~bm).sum())"
    r"\1PAIRHERE = PAIRW and abs(b - A.pairband) < 1e-9 \\"
    r"\1    and abs(frac - A.pairfrac) < 1e-9"
    r"\1if PAIRHERE and (q, b) not in PDONE:"
    r"\1    PDONE.add((q, b))"
    r"\1    res.append({'meth': 'pairdenom', 'ds': A.ds, 'q': q, 'b': b,"
    r"\1                'frac': -1.0, 'seed': -1, 'MI': MI, 'MO': MO,"
    r"\1                'pnT': pn(tm), 'pnB': pn(bm)})\n")

sub(r"'seed': sdd, 'M': int\(ar\.size\), 'MT': MT,",
    "'seed': sdd, 'M': int(ar.size), 'MT': MT, 'MI': MI, 'MO': MO,")

# ----------------------------------------------------- per-pair error counts --
sub(r"(\n                    )if extra: row\.update\(extra\)\n",
    r"\1if PAIRHERE and mt in PAIRM:"
    r"\1    row['pe'] = pv(pred, Ht, tm)"
    r"\1    row['qe'] = pv(pred, Ht, bm)"
    r"\1if extra: row.update(extra)\n")

# -------------------------------------------------------------- new designs --
sub(r"rec\('ig', post\(Lh, sg, cg, SIGSC, sd\) >= tau\)",
    "rec('ig', post(Lh, sg, cg, SIGSC, sd) >= tau)\n"
    "                rs2 = np.random.default_rng(5501 + 7919 * sdd\n"
    "                                            + int(1e5 * frac) + int(1e3 * q)\n"
    "                                            + int(1e2 * b))\n"
    "                if DESIGNS is None or 'ucb' in DESIGNS:\n"
    "                    su2, cu2 = greedy(rs2, tau, Lh, vhat, T, c1, s1, A.igrounds, 'ucb')\n"
    "                    rec('ucb', post(Lh, su2, cu2, SIGSC, sd) >= tau)\n"
    "                if DESIGNS is None or 'straddle' in DESIGNS:\n"
    "                    ss2, cs2 = greedy(rs2, tau, Lh, vhat, T, c1, s1, A.igrounds, 'straddle')\n"
    "                    rec('straddle', post(Lh, ss2, cs2, SIGSC, sd) >= tau)")

# ------------------------------------------------------------- output name ----
sub(r"OUTF = 'zz_", "OUTF = 'ww_")

open(DST, 'w').write(s)
if fails:
    print('FAILED:', fails)
    sys.exit(1)
import ast
ast.parse(s)
print('wrote %s, %d lines, parses clean' % (DST, s.count('\n') + 1))
