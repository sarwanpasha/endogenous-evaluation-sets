# Paired bootstrap over seeds (4000 resamples) of the inflation factor
#   x = mean_s[ R(UNIF,B) - R(D,B) ]  /  mean_s[ R(UNIF,B*) - R(D,B*) ]
# computed from the per-seed records in xy_{screen}_svd_homo.json at q = 0.90,
# T = 0.32N.  "n/a" means the design has no positive gain to inflate on one or
# both sets, so no ratio is defined.  20 seeds on O'Neil, 10 on DrugComb.
# gs = gain on B*, gb = gain on B.
INFLATION = {
 ('oneil', 0.30): {
   'meas':     dict(gs=-0.0076, gb=+0.0008, x=None, ci=None),
   'passive':  dict(gs=-0.0244, gb=-0.0424, x=None, ci=None),
   'band_hat': dict(gs=+0.0075, gb=+0.0335, x=4.44, ci=(2.92, 9.99)),
   'lse':      dict(gs=+0.0165, gb=+0.0448, x=2.72, ci=(2.17, 3.85)),
   'band_inf': dict(gs=+0.0780, gb=+0.2338, x=3.00, ci=(2.89, 3.12)),
 },
 ('oneil', 0.10): {
   'meas':     dict(gs=+0.0003, gb=+0.0043, x=None, ci=None),  # unstable, 55% of resamples
   'passive':  dict(gs=+0.0021, gb=-0.0348, x=None, ci=None),
   'band_hat': dict(gs=+0.0101, gb=-0.0120, x=None, ci=None),
   'lse':      dict(gs=+0.0224, gb=-0.0049, x=None, ci=None),
   'band_inf': dict(gs=+0.0307, gb=+0.3184, x=10.38, ci=(8.86, 12.51)),
 },
 ('dcomb', 0.30): {
   'meas':     dict(gs=-0.0174, gb=-0.0128, x=None, ci=None),
   'passive':  dict(gs=-0.0984, gb=-0.1175, x=None, ci=None),
   'band_hat': dict(gs=-0.0344, gb=-0.0299, x=None, ci=None),
   'lse':      dict(gs=-0.0349, gb=-0.0302, x=None, ci=None),
   'band_inf': dict(gs=+0.0242, gb=+0.0425, x=1.76, ci=(1.72, 1.81)),
 },
 ('dcomb', 0.10): {
   'meas':     dict(gs=-0.0082, gb=-0.0057, x=None, ci=None),
   'passive':  dict(gs=-0.0415, gb=-0.0771, x=None, ci=None),
   'band_hat': dict(gs=-0.0070, gb=-0.0142, x=None, ci=None),
   'lse':      dict(gs=-0.0034, gb=-0.0081, x=None, ci=None),
   'band_inf': dict(gs=+0.0155, gb=+0.1606, x=10.38, ci=(9.63, 11.49)),
 },
}
