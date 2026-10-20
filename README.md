# Endogenous Evaluation Sets in Adaptive Drug-Combination Screening

Code for the paper of the same name.

A combination screen spends a fixed budget of wells and returns a shortlist of
drug pairs whose synergy score exceeds a threshold. The decision is contested
only for pairs whose score lies near that threshold, so adaptive designs
concentrate measurement there. Reporting that this helps requires committing to
a set of pairs that count as "near the threshold", and the only quantity
available to locate that set before the budget is spent is a fitted model of the
score surface, which is the same model the design aims with. The evaluation set
is then endogenous to the design. This repository contains everything used to
measure how far that moves the reported numbers.

Two boundary sets appear throughout:

* `B*` (written `berrT` in the code) is defined by the withheld values,
  `{e : |F_e - tau| <= b * sigma_F}`. It is model free.
* `B` (written `berr` in the code) is defined by a completion `L`,
  `{e : |L_e - tau| <= b * sigma_F}`. It is what a practitioner can actually
  construct, because `F` is exactly what the screen does not know.

The inflation factor is a design's improvement over the non-targeting control on
`B`, divided by its improvement on the same control on `B*`.

## Layout

Everything sits at the top level, because `xy.py` and `x12.py` pull their shared
preludes in with `exec(open(...))` rather than `import` and so expect their
dependencies in the working directory.

### Replay experiments on the two public screens

| File | Role |
| --- | --- |
| `x1.py` | Shared prelude: loads a screen pickle into a pair-by-cell-line matrix (`BLD`) and fits the low-rank completion (`fit`). |
| `x12.py` | Shared prelude: budget allocation (`alloc`), the two-stage straddle allocation (`apt2`), and the fused posterior mean (`post`). |
| `feat.py` | Compound and cell-line features for the learned priors. |
| `xy.py` | The main driver. Replays the full design lattice against a fixed empirical ground truth and writes one record per (budget, seed, quantile, band, design). |
| `zy.py` | The extended driver. Adds the three measurements described below. |
| `mkpatch.py` | Builds `zy.py` from `xy.py`, so the difference between the two is auditable as a patch rather than a fork. |
| `wy.py` | The second extended driver. Adds the subset decomposition, the per-pair error counts and two further designs. |
| `mkpatch2.py` | Builds `wy.py` from `zy.py`, on the same principle. |

Run from the directory holding the screen pickles:

```
python3 xy.py --ds oneil --prior svd --noise homo --seeds 20 \
    --fracs 0.0025,0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12 \
    --bands 0.05,0.10,0.20,0.30,0.50 --quants 0.80,0.85,0.90,0.95,0.99
```

Results are appended to `xy_<screen>_<score>_<prior>_<noise>.json` and the run
resumes from whatever is already there, so a job that hits a wall-clock limit can
simply be resubmitted.

Useful flags: `--prior` selects the prior class (`svd` for the low-rank
completion, `gbt` for the boosted regressor, `clf` and `clfu` for the classifier
trained on the sign at the threshold); `--noise` selects `homo` for the
optimistic homoskedastic model or `emp` for the model calibrated to each
screen's own repeated wells; `--a1` sets the two-stage split; `--rank` overrides
the completion rank; `--designs` restricts which arms are run.

### What `zy.py` adds

`zy.py` is `xy.py` plus three measurements. Regenerate it with
`python3 mkpatch.py`, which asserts that every edit site matched, so the patch
fails loudly rather than silently if the driver changes.

1. **`berrH`** is the boundary error on a band built from `Lh`, the stage-one
   prior fit on purchased wells only, instead of from `Lv`, the reference
   completion. `Lv` is a truncated SVD of the complete matrix and so has seen
   the values the designs are scored against; it bounds the effect but is not a
   band anyone could construct. `berrH` is the realistic counterpart.
2. **`hit100`, `hit500`, `hitpos`** are the precision of the top-k wells ranked
   by the fused score, with `hitpos` evaluated at k equal to the number of true
   positives. This is the hit-discovery endpoint that deployed adaptive screens
   report, as opposed to a zero-one rate on a band.
3. **Cost-aware arms** (`*_k3`, `*_k10`) threshold the posterior probability at
   `1/(1+k)` instead of thresholding the fused mean at `tau`. Under an
   asymmetric cost that prices a missed combination at `k` wasted confirmatory
   assays, this is the Bayes rule, so the cost comparison is no longer being
   made under a rule that ignores the cost.

The random seed stream is untouched, so the arms that exist in both drivers
reproduce identical numbers.

### What `wy.py` adds

`wy.py` is `zy.py` plus four things. Regenerate it with `python3 mkpatch2.py`,
which asserts that every edit site matched.

1. **Subset error rates.** Every arm is scored separately on the overlap
   `B n B*`, on the part of the model-free band that the model-defined band
   misses, `B* \ B`, and on `B` itself. Against the non-targeting control these
   give the overlap `omega` and the two corrections `kappa` and `mu` that the
   ratio proposition names, so the decomposition is measured rather than
   assumed.
2. **Per-pair error counts.** Errors on `B*` and on `B` are summed within each
   drug pair and written as `pe` and `qe`, with the matching denominators
   written once per band as a `pairdenom` record. An interval can then be
   resampled over drug pairs, which is the unit a screen can add or drop, rather
   than over noise seeds.
3. **`ucb`**, which spends the second stage where the posterior mean plus two
   standard deviations is largest. This is the exploitation loop that deployed
   recommenders run.
4. **`straddle`**, which spends the second stage where `1.96 s(e) - |m(e) -
   tau|` is largest. This is the level-set acquisition of the active level-set
   estimation literature, applied to the same posterior.

The two new designs draw from their own generator, so every arm that also
exists in `zy.py` reproduces identical numbers.

```
python3 wy.py --ds oneil --prior svd --noise homo --seeds 20 \
    --fracs 0.08,0.32,1.28 --bands 0.10,0.30 --quants 0.90 --tag _p1
```

Records are appended to `ww_<screen>_<score>_<prior>_<noise><tag>.json`. The
per-pair vectors are long, so run `wy.py` on the settings the paper reports
rather than on the full lattice; `--pairs 0` turns them off.

### The controlled study

`syn.py` builds the situation directly, with no screen involved. A rank-`r`
matrix is drawn and standardised, `tau` is its `q`-quantile, and a reference
prior `L = F + eta * xi` fixes the model-defined band. A design spends a fraction
`alpha` of its budget uniformly inside that band and the rest uniformly over the
support, so `alpha = 0` is the non-targeting control and `alpha = 1` is the
design whose allocation is the band itself.

```
python3 syn.py --P 800 --m 80 --rank 8 --eta 0.50 --sig 0.50 \
    --q 0.90 --b 0.30 --frac 0.32 --seeds 120 \
    --alphas 0,0.25,0.5,0.75,1.0 --sweep base --out syn_base.json
```

Repeat with `--sweep b`, `--sweep eta`, `--sweep rank` and `--sweep frac`, then
concatenate the five JSON files into one list. The estimation noise is not a
confound here: `lambda` is set to its exact value `sig**2 / eta**2` rather than
estimated, so what the sweeps move is the targeting and nothing else.

Defaults used in the paper: `P = 800`, `m = 80`, `r = 8`, `eta = 0.50 sigma_F`,
`sigma = 0.50 sigma_F`, `q = 0.90`, `b = 0.30`, `T = 0.32 N`, 120 seeds per
point. Standard errors over seeds are at most 0.002 on either boundary set.

### Summarisation, intervals and figures

| File | Role |
| --- | --- |
| `sy.py` | Summarises every `xy_*.json` into CSVs, with bootstrap intervals over seeds on boundary error and global F1. |
| `slim.py` | Cuts the full summary down to the columns the paper uses. |
| `mko.py` | Turns the summary CSVs into the numbers and tables quoted in the paper. |
| `bootstrap_ci.py` | Paired bootstrap intervals on the inflation factor, over seeds. `--glob` selects which record files to read. |
| `pair_ci.py` | Intervals on the two gains and on their difference, resampled over drug pairs, with the measured `omega`, `kappa` and `mu` beside the predicted factor. |
| `alpha_decomp.py` | The same decomposition in the controlled model, swept over targeting strength. |
| `band_compare.py` | The three-way band comparison, the hit rates and the cost-aware arms, from the `zz_*.json` files written by `zy.py`. |
| `figs.py` | Renders the three figures from `series.py`. |

`bootstrap_ci.py` resamples seeds and recomputes the whole ratio, because the
factor is a ratio of two averages and an interval cannot be assembled from
separate intervals on the numerator and the denominator. Resamples that leave a
non-positive gain on either set yield no ratio; if fewer than 90 percent of
resamples are usable the factor is reported as undefined rather than quoted from
the usable tail.

```
python3 sy.py                       # run where the xy_*.json files are
python3 bootstrap_ci.py --dir . --frac 0.32 --q 0.90
python3 band_compare.py zz_oneil_Y_svd_homo_v2.json --frac 0.32
python3 pair_ci.py ww_oneil_Y_svd_homo_p1.json --frac 0.32 --q 0.90 --b 0.30
python3 alpha_decomp.py --seeds 200 --out alpha_base.json
python3 figs.py                     # writes figures/{bands,effect,budget}.pdf
                                    # set FIGDIR to write them elsewhere
```

### Recorded series

`series.py` holds the exact series plotted in the paper, with the provenance of
each block noted inline, so the figures can be regenerated without re-running
any experiment. `ci.py` holds the bootstrap intervals in the same form.
`hat.py` holds the three-way band comparison, the hit rates and the
cost-aware costs for both screens. `alpha.py` holds the measured decomposition
in the controlled model, `robust.py` the factors under the four synergy models,
and `pairs.py` the gains and their difference with intervals over drug pairs. `test_small.npz` is a synthetic 60-by-40 screen with 85 percent of its cells
observed, for exercising the pipeline end to end without downloading anything.

## The screens

Both screens are public and are obtained through the Therapeutics Data Commons.

| Screen | Pairs | Cell lines | Measured wells | Score | Completion rank |
| --- | --- | --- | --- | --- | --- |
| O'Neil | 583 | 39 | 22,737 | `Y` | 8 |
| DrugComb | 5,427 | 59 | 287,268 | `Synergy_Loewe` | 12 |

The drivers expect `df_OncoPolyPharmacology.pkl` and `df_DrugComb.pkl` in the
working directory: pandas DataFrames with columns `Drug1_ID`, `Drug2_ID`,
`Cell_Line_ID` and the score column. `x1.py` pivots these into a pair-by-cell-line
matrix, averaging wells measured more than once, which is also how the replicate
noise model is calibrated.

## Requirements

Python 3.9 or later, with `numpy`, `scipy`, `pandas`, `scikit-learn` and
`matplotlib`. The learned priors use scikit-learn's gradient boosting; nothing
needs a GPU.

```
pip install numpy scipy pandas scikit-learn matplotlib
```

## Reproducing the paper

1. Place the two screen pickles beside the scripts.
2. Run `xy.py` for each screen under each prior class and noise model.
3. Run `zy.py` for each screen to add the purchased-data band, the hit rates and
   the cost-aware arms, and `wy.py` for the subset decomposition, the per-pair
   counts and the two further designs.
4. Run `sy.py`, then `slim.py`, then `mko.py`.
5. Run `syn.py` over the five sweeps.
6. Run `bootstrap_ci.py` and `pair_ci.py` for the intervals, `alpha_decomp.py`
   for the decomposition, and `figs.py` for the figures.

Every design is replayed against the same fixed empirical ground truth with a
common random seed stream, so comparisons at a fixed budget are paired across
designs.
