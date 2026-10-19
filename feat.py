import numpy as np, pandas as pd

def load(p, scol):
 d = pd.read_pickle(p)
 a = d['Drug1_ID'].astype(str).values; b = d['Drug2_ID'].astype(str).values
 k = np.where(a < b, a + '|' + b, b + '|' + a)
 c = d['Cell_Line_ID'].astype(str).values
 y = pd.to_numeric(d[scol], errors='coerce').values
 pv = pd.DataFrame({'k': k, 'c': c, 'y': y}).pivot_table(index='k', columns='c', values='y', aggfunc='mean')
 F = pv.values.astype(float)
 keys = np.asarray(pv.index); cells = np.asarray(pv.columns)
 sm = {}
 if 'Drug1' in d.columns:
  z1 = pd.DataFrame({'i': a, 's': d['Drug1'].astype(str).values}).drop_duplicates('i')
  z2 = pd.DataFrame({'i': b, 's': d['Drug2'].astype(str).values}).drop_duplicates('i')
  sm = dict(zip(z1['i'], z1['s'])); sm.update(dict(zip(z2['i'], z2['s'])))
 ce = {}
 if 'CellLine' in d.columns:
  vals = d['CellLine'].values
  for i in range(len(c)):
   cid = c[i]
   if cid in ce: continue
   x = vals[i]
   if isinstance(x, (list, tuple)):
    v = np.concatenate([np.asarray(z).ravel().astype(np.float32) for z in x])
   else:
    v = np.asarray(x).ravel().astype(np.float32)
   ce[cid] = v
 del d
 return F, keys, cells, sm, ce
