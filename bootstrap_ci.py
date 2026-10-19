"""Paired bootstrap intervals for the inflation factor.

The inflation factor of a design D is

    x(D) = [ R(UNIF, B) - R(D, B) ] / [ R(UNIF, B*) - R(D, B*) ],

the improvement D makes over the non-targeting control on the model-defined
band, divided by the improvement it makes on the model-free band.  Both
numerator and denominator are averages over seeds, so the interval has to be
formed by resampling seeds and recomputing the ratio, not by combining separate
intervals on the two gains.  Seeds are resampled jointly for the two sets, which
keeps the pairing that makes the comparison meaningful.

A resample that produces a non-positive denominator or numerator yields no
ratio.  The fraction of usable resamples is reported; below 90 percent we treat
the factor as undefined rather than quoting an interval built from the usable
tail, because that tail is not a sample from the quantity of interest.

Usage:

    python3 bootstrap_ci.py --dir /path/to/results --frac 0.32 --q 0.90

Expects the per-seed record files written by screens/xy.py, named
xy_<screen>_<score>_<prior>_<noise>.json.
"""
import argparse
import glob
import json
import os

import numpy as np

REFERENCE = 'blind'            # the non-targeting control, UNIF in the paper
DESIGNS = ['meas', 'passive', 'band_hat', 'lse', 'ig', 'band_inf']
RESAMPLES = 4000
MIN_GAIN = 1e-9                # a gain at or below this is treated as absent
MIN_USABLE = 0.90              # fraction of resamples that must yield a ratio


def load(path, q):
    """Index the per-seed records by (band, budget, design, seed)."""
    table = {}
    for row in json.load(open(path)):
        if abs(float(row['q']) - q) > 1e-9:
            continue
        key = (float(row['b']), float(row['frac']), row['meth'], int(row['seed']))
        table[key] = (float(row['berrT']), float(row['berr']))
    return table


def gains(table, design, band, frac, seeds):
    """Per-seed improvement of `design` over the control, on B* and on B."""
    star, model = [], []
    for s in seeds:
        ref = table.get((band, frac, REFERENCE, s))
        cur = table.get((band, frac, design, s))
        if ref is None or cur is None:
            return None, None
        star.append(ref[0] - cur[0])
        model.append(ref[1] - cur[1])
    return np.asarray(star), np.asarray(model)


def interval(star, model, rng):
    point = (model.mean() / star.mean()
             if star.mean() > MIN_GAIN and model.mean() > 0 else None)
    draws = []
    n = star.size
    for _ in range(RESAMPLES):
        idx = rng.integers(0, n, n)
        a, b = star[idx].mean(), model[idx].mean()
        if a > MIN_GAIN and b > 0:
            draws.append(b / a)
    usable = len(draws) / RESAMPLES
    if usable < MIN_USABLE:
        return point, None, usable
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return point, (float(lo), float(hi)), usable


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='.', help='directory holding xy_*.json')
    ap.add_argument('--q', type=float, default=0.90)
    ap.add_argument('--frac', type=float, default=0.32)
    ap.add_argument('--bands', default='0.30,0.10')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    files = sorted(glob.glob(os.path.join(args.dir, 'xy_*_svd_homo.json')))
    if not files:
        raise SystemExit('no xy_*_svd_homo.json under %s' % args.dir)

    for path in files:
        table = load(path, args.q)
        seeds = sorted({k[3] for k in table})
        print('== %s  (%d seeds)' % (os.path.basename(path), len(seeds)))
        for band in [float(x) for x in args.bands.split(',')]:
            for design in DESIGNS:
                star, model = gains(table, design, band, args.frac, seeds)
                if star is None:
                    continue
                point, ci, usable = interval(star, model, rng)
                shown = '%.2f' % point if point is not None else 'none'
                span = '[%.2f, %.2f]' % ci if ci else 'none'
                print('   b=%.2f %-9s gain* %+.4f  gain %+.4f  factor %-6s %s'
                      '  (usable %.0f%%)'
                      % (band, design, star.mean(), model.mean(), shown, span,
                         100 * usable))


if __name__ == '__main__':
    main()
