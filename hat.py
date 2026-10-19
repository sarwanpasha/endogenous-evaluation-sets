# Results from zz_oneil_Y_svd_homo_v2.json (driver zy.py), O'Neil, prior svd,
# homoskedastic noise, q = 0.90 only, 20 seeds, T in {0.04, 0.32, 1.28} N.
#
# This is a separate run from the one behind Table 1: it sweeps a single
# quantile, so the random stream differs and the numbers are not identical to
# the five-quantile run.  Everything quoted together below comes from this one
# run, so the comparison inside it is internally consistent.
#
# errT = error on B* (model free), errB = error on B (reference completion Lv,
# fit to the complete matrix), errH = error on Bhat (stage-one fit Lh, fit on
# purchased wells only).  xB and xH are the inflation factors against B and
# against Bhat, with paired bootstrap intervals over seeds (4000 resamples).

BANDS = {}  # band -> design -> dict
BANDS[0.05] = dict(MT=340, MB=287, MH=304, arms={
 'meas':     dict(errT=0.5066, errB=0.3974, errH=0.1868, xB=1.88, ciB=(0.67,13.66), xH=None, ciH=None),
 'passive':  dict(errT=0.4984, errB=0.4294, errH=0.2291, xB=None, ciB=None, xH=None, ciH=None),
 'band_hat': dict(errT=0.5015, errB=0.4286, errH=0.0403, xB=None, ciB=None, xH=14.84, ciH=(6.74,121.24)),
 'lse':      dict(errT=0.4972, errB=0.4261, errH=0.0203, xB=None, ciB=None, xH=11.68, ciH=(6.46,52.43)),
 'ig':       dict(errT=0.5063, errB=0.4472, errH=0.3467, xB=None, ciB=None, xH=None, ciH=None),
 'band_inf': dict(errT=0.4554, errB=0.0617, errH=0.1242, xB=6.20, ciB=(5.43,7.14), xH=1.01, ciH=(0.76,1.30)),
})
BANDS[0.10] = dict(MT=678, MB=583, MH=652, arms={
 'meas':     dict(errT=0.4779, errB=0.3961, errH=0.1842, xB=14.54, ciB=None, xH=None, ciH=None),
 'passive':  dict(errT=0.4761, errB=0.4352, errH=0.2265, xB=None, ciB=None, xH=None, ciH=None),
 'band_hat': dict(errT=0.4669, errB=0.4123, errH=0.0633, xB=None, ciB=None, xH=10.28, ciH=(6.23,36.10)),
 'lse':      dict(errT=0.4577, errB=0.4050, errH=0.0325, xB=None, ciB=None, xH=7.18, ciH=(5.28,12.21)),
 'ig':       dict(errT=0.4745, errB=0.4454, errH=0.3887, xB=None, ciB=None, xH=None, ciH=None),
 'band_inf': dict(errT=0.4454, errB=0.0804, errH=0.1118, xB=9.77, ciB=(8.48,11.61), xH=2.06, ciH=(1.69,2.65)),
})
BANDS[0.20] = dict(MT=1394, MB=1222, MH=1332, arms={
 'meas':     dict(errT=0.4221, errB=0.3936, errH=0.1802, xB=None, ciB=None, xH=None, ciH=None),
 'passive':  dict(errT=0.4290, errB=0.4270, errH=0.2219, xB=None, ciB=None, xH=None, ciH=None),
 'band_hat': dict(errT=0.4095, errB=0.3727, errH=0.0865, xB=2.83, ciB=(1.75,6.67), xH=10.60, ciH=(6.84,23.88)),
 'lse':      dict(errT=0.3925, errB=0.3596, errH=0.0564, xB=1.45, ciB=(1.18,1.87), xH=4.70, ciH=(3.85,6.13)),
 'ig':       dict(errT=0.4251, errB=0.4364, errH=0.3826, xB=None, ciB=None, xH=None, ciH=None),
 'band_inf': dict(errT=0.3548, errB=0.1141, errH=0.0923, xB=4.47, ciB=(4.22,4.75), xH=1.32, ciH=(1.20,1.47)),
})
BANDS[0.30] = dict(MT=2164, MB=1870, MH=2048, arms={
 'meas':     dict(errT=0.3802, errB=0.3683, errH=0.1741, xB=None, ciB=None, xH=None, ciH=None),
 'passive':  dict(errT=0.3971, errB=0.4115, errH=0.2174, xB=None, ciB=None, xH=None, ciH=None),
 'band_hat': dict(errT=0.3635, errB=0.3348, errH=0.1122, xB=3.74, ciB=(2.70,7.02), xH=6.29, ciH=(4.36,11.70)),
 'lse':      dict(errT=0.3577, errB=0.3247, errH=0.0937, xB=2.97, ciB=(2.40,4.03), xH=5.10, ciH=(4.15,6.92)),
 'ig':       dict(errT=0.3858, errB=0.4107, errH=0.3564, xB=None, ciB=None, xH=None, ciH=None),
 'band_inf': dict(errT=0.2948, errB=0.1357, errH=0.0833, xB=3.00, ciB=(2.90,3.09), xH=1.11, ciH=(1.05,1.19)),
})
BANDS[0.50] = dict(MT=3849, MB=3488, MH=3768, arms={
 'meas':     dict(errT=0.2976, errB=0.3119, errH=0.1625, xB=None, ciB=None, xH=None, ciH=None),
 'passive':  dict(errT=0.3319, errB=0.3661, errH=0.2065, xB=None, ciB=None, xH=None, ciH=None),
 'band_hat': dict(errT=0.2966, errB=0.3060, errH=0.2117, xB=None, ciB=None, xH=None, ciH=None),
 'lse':      dict(errT=0.2947, errB=0.3040, errH=0.2119, xB=None, ciB=None, xH=None, ciH=None),
 'ig':       dict(errT=0.3122, errB=0.3537, errH=0.2951, xB=None, ciB=None, xH=None, ciH=None),
 'band_inf': dict(errT=0.2097, errB=0.1430, errH=0.0730, xB=2.15, ciB=(2.10,2.20), xH=1.11, ciH=(1.07,1.16)),
})

# Precision of the top-k wells ranked by the fused score, b = 0.30.
HITS = {
 0.04: {'blind':(0.6720,0.1835,0.0669), 'passive':(0.5765,0.3553,0.2073),
        'band_hat':(0.7840,0.3158,0.1657), 'lse':(0.7890,0.3162,0.1658),
        'band_inf':(1.0000,0.9560,0.6822)},
 0.32: {'blind':(0.9995,0.8023,0.2712), 'passive':(0.5450,0.4980,0.3854),
        'band_hat':(0.8365,0.7249,0.4028), 'lse':(0.8345,0.7205,0.4020),
        'band_inf':(1.0000,0.9750,0.7235)},
 1.28: {'blind':(1.0000,0.9932,0.6247), 'passive':(0.6725,0.7522,0.6619),
        'band_hat':(0.6040,0.6050,0.6338), 'lse':(0.5795,0.5827,0.6595),
        'band_inf':(1.0000,0.9771,0.7290)},
}

# Asymmetric cost on B* (b = 0.30) under the cost-blind point rule and under
# the Bayes rule that thresholds the posterior at 1/(1+k).
COSTS = {
 0.04: {'blind':((1.133,1.125),(3.767,0.631)), 'passive':((1.103,1.103),(3.623,3.615)),
        'band_hat':((1.121,1.107),(3.709,0.897)), 'band_inf':((0.825,0.663),(2.539,0.893))},
 0.32: {'blind':((1.049,1.007),(3.415,0.682)), 'passive':((0.900,0.898),(2.659,2.638)),
        'band_hat':((0.937,0.887),(2.943,1.228)), 'band_inf':((0.715,0.626),(2.185,0.884))},
 1.28: {'blind':((0.833,0.731),(2.564,0.747)), 'passive':((0.767,0.649),(2.271,1.159)),
        'band_hat':((0.710,0.601),(2.074,1.325)), 'band_inf':((0.633,0.566),(1.961,0.827))},
}

# ---------------------------------------------------------------- DrugComb --
# zz_dcomb_Synergy_Loewe_svd_homo_v2.json, same configuration, 10 seeds.
# No feasible targeting design has a positive gain on either set at any band,
# so only ORACLE admits a factor.
DC_BANDS = {
 0.05: dict(MT=12298, MB=7708,   MH=7165,
            blind=(0.4479, 0.3730, 0.1401),
            oracle=(0.4438, 0.1442, 0.1081),
            gains=(0.0041, 0.2288, 0.0320),
            xB=56.17, ciB=(43.43, 81.21), xH=7.86, ciH=(6.05, 11.43)),
 0.10: dict(MT=25094, MB=15878,  MH=14325,
            blind=(0.4017, 0.3531, 0.1407),
            oracle=(0.3850, 0.1907, 0.1036),
            gains=(0.0167, 0.1624, 0.0371),
            xB=9.71, ciB=(9.24, 10.24), xH=2.22, ciH=(2.07, 2.38)),
 0.20: dict(MT=53562, MB=36314,  MH=29136,
            blind=(0.3176, 0.3088, 0.1393),
            oracle=(0.2916, 0.2271, 0.1020),
            gains=(0.0259, 0.0817, 0.0373),
            xB=3.15, ciB=(3.07, 3.26), xH=1.44, ciH=(1.37, 1.52)),
 0.30: dict(MT=83866, MB=63300,  MH=45147,
            blind=(0.2578, 0.2568, 0.1377),
            oracle=(0.2334, 0.2139, 0.1058),
            gains=(0.0244, 0.0429, 0.0320),
            xB=1.76, ciB=(1.73, 1.80), xH=1.31, ciH=(1.27, 1.36)),
 0.50: dict(MT=136613, MB=121159, MH=81511,
            blind=(0.1944, 0.1924, 0.1327),
            oracle=(0.1703, 0.1653, 0.1073),
            gains=(0.0241, 0.0271, 0.0254),
            xB=1.13, ciB=(1.12, 1.14), xH=1.06, ciH=(1.01, 1.10)),
}

DC_HITS = {
 0.04: {'blind':(0.8590,0.6058,0.1057), 'passive':(0.3740,0.3434,0.1500),
        'band_hat':(0.5350,0.4900,0.1370), 'lse':(0.5380,0.4826,0.1370),
        'band_inf':(1.0000,0.9640,0.4725)},
 0.32: {'blind':(0.9970,0.9054,0.2391), 'passive':(0.2340,0.2334,0.2275),
        'band_hat':(0.3840,0.3242,0.2735), 'lse':(0.3840,0.3260,0.2743),
        'band_inf':(1.0000,0.9640,0.5356)},
 1.28: {'blind':(1.0000,0.9910,0.4402), 'passive':(0.6370,0.5902,0.4482),
        'band_hat':(0.2950,0.2736,0.4220), 'lse':(0.2720,0.2676,0.4114),
        'band_inf':(1.0000,0.9648,0.6110)},
}

DC_COSTS = {
 0.04: {'blind':((0.713,0.713),(2.362,0.763)), 'passive':((0.709,0.709),(2.288,2.285)),
        'band_hat':((0.710,0.711),(2.338,0.797)), 'band_inf':((0.610,0.563),(1.905,1.030))},
 0.32: {'blind':((0.678,0.679),(2.150,0.771)), 'passive':((0.670,0.670),(1.767,1.754)),
        'band_hat':((0.661,0.674),(1.952,1.019)), 'band_inf':((0.555,0.526),(1.681,1.028))},
 1.28: {'blind':((0.589,0.595),(1.674,0.766)), 'passive':((0.568,0.577),(1.520,0.889)),
        'band_hat':((0.553,0.536),(1.432,1.024)), 'band_inf':((0.477,0.462),(1.439,0.978))},
}
