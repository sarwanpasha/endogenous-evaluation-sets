"""Summarise all xy_*.json result files into one compact CSV for the paper."""
import json, glob, os, csv, sys
import numpy as np

B = 600                      # bootstrap resamples
CI = ('berrT', 'f1', 'cost3', 'berrBal', 'gcost3')
MEAN = ('berr', 'gerr', 'prec', 'rec', 'bprec', 'brec', 'bfp', 'bfn', 'cost10',
        'gfp', 'gfn', 'gcost10', 'nsel', 'posfrac', 'pri_e', 'prip_e', 'prih_e',
        'M', 'MT', 'Mhat', 'lam', 'pred')
KEYS = ('ds', 'score', 'prior', 'noise', 'q', 'b', 'frac', 'meth')

def boot(v, rng):
    v = np.asarray(v, float)
    if v.size < 2:
        return float(v.mean()), float('nan'), float('nan')
    idx = rng.integers(0, v.size, (B, v.size))
    m = v[idx].mean(1)
    return float(v.mean()), float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))

rows = []
files = sorted(glob.glob('xy_*.json'))
meta = []
for f in files:
    try:
        d = json.load(open(f))
    except Exception as e:
        print('SKIP', f, e, file=sys.stderr); continue
    src = os.path.basename(f)[:-5]
    for r in d:
        r['srcname'] = src
        rows.append(r)
    meta.append((src, len(d)))
print(json.dumps({'files': meta, 'rows': len(rows)}), flush=True)

grp = {}
for r in rows:
    k = (r['srcname'],) + tuple(r.get(x) for x in KEYS)
    grp.setdefault(k, []).append(r)

rng = np.random.default_rng(7)
out = []
for k, v in sorted(grp.items(), key=lambda kv: [str(x) for x in kv[0]]):
    o = dict(zip(('srcname',) + KEYS, k))
    o['n'] = len(v)
    for mt in CI:
        vals = [x[mt] for x in v if x.get(mt) is not None]
        if vals:
            mu, lo, hi = boot(vals, rng)
            o[mt] = round(mu, 5); o[mt + '_lo'] = round(lo, 5); o[mt + '_hi'] = round(hi, 5)
    for mt in MEAN:
        vals = [x[mt] for x in v if x.get(mt) is not None]
        if vals:
            o[mt] = round(float(np.mean(vals)), 5)
    out.append(o)

cols = ['srcname'] + list(KEYS) + ['n']
for mt in CI:
    cols += [mt, mt + '_lo', mt + '_hi']
cols += list(MEAN)
w = csv.DictWriter(open('summary_y.csv', 'w', newline=''), fieldnames=cols, extrasaction='ignore')
w.writeheader()
for o in out:
    w.writerow(o)
print(json.dumps({'k': 'wrote', 'file': 'summary_y.csv', 'groups': len(out)}), flush=True)

# paper-scoped slim version: default q, the three headline budgets, all designs,
# plus the full frac grid at the default band
def keep(o):
    q = float(o['q']); b = float(o['b']); fr = float(o['frac'])
    if abs(q - 0.90) < 1e-9 and b in (0.10, 0.30):
        return True
    if abs(b - 0.30) < 1e-9:
        return True
    if fr in (0.04, 0.32, 1.28):
        return True
    return False

slim = [o for o in out if keep(o)]
w = csv.DictWriter(open('mini_y.csv', 'w', newline=''), fieldnames=cols, extrasaction='ignore')
w.writeheader()
for o in slim:
    w.writerow(o)
print(json.dumps({'k': 'wrote', 'file': 'mini_y.csv', 'groups': len(slim)}), flush=True)
print('SUMMARY_DONE', flush=True)
