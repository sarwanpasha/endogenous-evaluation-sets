"""Compare the three boundary sets, the hit rates and the cost-aware arms.

Reads the per-seed records written by zy.py (zz_*.json) and reports, for each
band half-width at one budget:

  * the error on B*   (model free, from the withheld values),
  * the error on B    (from the reference completion, fit to the full matrix),
  * the error on Bhat (from the stage-one prior, fit on purchased wells only),
  * the inflation factor against B and against Bhat, each with a paired
    bootstrap interval over seeds.

It then reports the precision of the top-k wells ranked by the fused score, and
the asymmetric cost under the cost-blind point rule against the Bayes rule for
that cost.

The point of the three-way comparison is that a design is flattered by the band
belonging to the model it aims with: the oracle design, which aims with the
reference completion, is inflated on B and not on Bhat, while the two-stage
designs, which aim with the stage-one prior, are inflated on Bhat.

Usage:

    python3 band_compare.py zz_oneil_Y_svd_homo_v2.json --frac 0.32
"""
import argparse
import json

import numpy as np

REFERENCE = 'blind'
DESIGNS = ['meas', 'passive', 'band_hat', 'lse', 'ig', 'band_inf']
HIT_DESIGNS = ['blind', 'passive', 'band_hat', 'lse', 'band_inf']
COST_DESIGNS = ['blind', 'passive', 'band_hat', 'band_inf']
RESAMPLES = 4000
MIN_GAIN = 1e-9
MIN_USABLE = 0.90


def index(path):
    table = {}
    for row in json.load(open(path)):
        table[(float(row['b']), float(row['frac']), row['meth'],
               int(row['seed']))] = row
    return table


def column(table, band, frac, design, seeds, key):
    vals = [table[(band, frac, design, s)].get(key) for s in seeds]
    vals = [v for v in vals if v is not None]
    return np.asarray(vals, dtype=float) if vals else None


def factor(star, other, rng):
    point = (other.mean() / star.mean()
             if star.mean() > MIN_GAIN and other.mean() > 0 else None)
    draws = []
    for _ in range(RESAMPLES):
        idx = rng.integers(0, star.size, star.size)
        a, b = star[idx].mean(), other[idx].mean()
        if a > MIN_GAIN and b > 0:
            draws.append(b / a)
    if len(draws) < MIN_USABLE * RESAMPLES:
        return point, None
    return point, tuple(np.round(np.quantile(draws, [0.025, 0.975]), 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('path')
    ap.add_argument('--frac', type=float, default=0.32)
    ap.add_argument('--cost-band', type=float, default=0.30)
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    table = index(args.path)
    seeds = sorted({k[3] for k in table})
    bands = sorted({k[0] for k in table})
    print('%s: %d seeds' % (args.path, len(seeds)))

    for band in bands:
        head = table[(band, args.frac, REFERENCE, seeds[0])]
        print('-- b=%.2f  |B*|=%d  |B|=%d  |Bhat|=%d'
              % (band, head['MT'], head['M'], head.get('MH', 0)))
        for design in DESIGNS:
            if (band, args.frac, design, seeds[0]) not in table:
                continue
            ref = [column(table, band, args.frac, REFERENCE, seeds, k)
                   for k in ('berrT', 'berr', 'berrH')]
            cur = [column(table, band, args.frac, design, seeds, k)
                   for k in ('berrT', 'berr', 'berrH')]
            if any(c is None for c in cur):
                continue
            star, on_b, on_h = (r - c for r, c in zip(ref, cur))
            xb, cb = factor(star, on_b, rng)
            xh, ch = factor(star, on_h, rng)
            print('   %-9s errT=%.4f errB=%.4f errH=%.4f | xB=%-6s %-14s'
                  ' | xH=%-6s %s'
                  % (design, cur[0].mean(), cur[1].mean(), cur[2].mean(),
                     ('%.2f' % xb) if xb else 'n/a', cb or '',
                     ('%.2f' % xh) if xh else 'n/a', ch or ''))

    print('\nprecision of the top k by fused score, b=%.2f' % args.cost_band)
    print('%-10s %6s %8s %8s %8s' % ('design', 'T/N', '@100', '@500', '@pos'))
    for frac in sorted({k[1] for k in table}):
        for design in HIT_DESIGNS:
            if (args.cost_band, frac, design, seeds[0]) not in table:
                continue
            vals = [column(table, args.cost_band, frac, design, seeds, k)
                    for k in ('hit100', 'hit500', 'hitpos')]
            if any(v is None for v in vals):
                continue
            print('%-10s %6.2f %8.4f %8.4f %8.4f'
                  % (design, frac, vals[0].mean(), vals[1].mean(),
                     vals[2].mean()))

    print('\ncost on B* under the point rule and the Bayes rule, b=%.2f'
          % args.cost_band)
    for frac in sorted({k[1] for k in table}):
        for design in COST_DESIGNS:
            parts = []
            for k in (3, 10):
                aware = '%s_k%d' % (design, k)
                if (args.cost_band, frac, aware, seeds[0]) not in table:
                    continue
                blind = column(table, args.cost_band, frac, design, seeds,
                               'cost%d' % k)
                smart = column(table, args.cost_band, frac, aware, seeds,
                               'cost%d' % k)
                parts.append('C%-2d point=%.3f bayes=%.3f'
                             % (k, blind.mean(), smart.mean()))
            if parts:
                print('T/N=%.2f %-10s %s' % (frac, design, '   '.join(parts)))


if __name__ == '__main__':
    main()
