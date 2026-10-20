"""Inflation factors on DrugComb under the four synergy models.

Each run replays the same design lattice against the same screen with the score
column changed, at q = 0.90 and T = 0.32N, ten seeds, paired bootstrap over
seeds.  Loewe is the model used in the body of the paper and its row is the run
behind the headline table; the other three come from the runs tagged _rb.

Entries are (factor, low, high) or None where the design has no gain on the
model-free band to inflate, so no ratio exists.
"""

ORACLE = {
 'Loewe': {0.05: (56.17, 43.0, 81.0), 0.10: (9.71, 9.2, 10.0),
           0.20: (3.15, 3.1, 3.3), 0.30: (1.76, 1.7, 1.8),
           0.50: (1.13, 1.1, 1.1)},
 'Bliss': {0.05: (147.52, 65.72, 1269.22), 0.10: (24.08, 21.65, 26.89),
           0.20: (7.67, 7.42, 7.93), 0.30: (4.88, 4.78, 4.98),
           0.50: (2.89, 2.87, 2.92)},
 'HSA':   {0.05: (47.06, 32.72, 95.38), 0.10: (20.29, 17.51, 24.01),
           0.20: (7.82, 7.57, 8.08), 0.30: (4.33, 4.28, 4.37),
           0.50: (2.46, 2.44, 2.48)},
 'ZIP':   {0.05: (59.74, 43.81, 94.49), 0.10: (23.41, 20.45, 27.58),
           0.20: (9.83, 9.20, 10.51), 0.30: (5.45, 5.31, 5.57),
           0.50: (3.02, 2.95, 3.08)},
}

ADAPT = {
 'Loewe': {0.05: None, 0.10: None, 0.20: None, 0.30: None, 0.50: None},
 'Bliss': {0.05: None, 0.10: None, 0.20: (1.07, 0.83, 1.37),
           0.30: (4.67, 3.74, 6.44), 0.50: None},
 'HSA':   {0.05: None, 0.10: None, 0.20: (1.18, 1.08, 1.27),
           0.30: (3.93, 3.42, 4.67), 0.50: None},
 'ZIP':   {0.05: None, 0.10: None, 0.20: (0.43, 0.30, 0.56),
           0.30: (3.22, 2.68, 4.08), 0.50: None},
}

# Gains on the model-free band for ADAPT at b = 0.30, which is what the factor
# above is multiplying; they are small, and the paper says so.
ADAPT_GAIN_STAR = {'Loewe': None, 'Bliss': 0.0065, 'HSA': 0.0088, 'ZIP': 0.0069}

BAND = {
 'HSA': {0.20: (4.91, 3.36, 17.78)},
}
