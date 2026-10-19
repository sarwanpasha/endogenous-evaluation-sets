import csv, json, os
KEEP = ['berrT','berrT_lo','berrT_hi','f1','f1_lo','f1_hi','berr','prec','pred','pri_e']
BASE = ['src','ds','score','prior','q','b','frac','meth','M','MT','n']
rows = list(csv.DictReader(open('summary.csv')))
srcs = sorted({r['src'] for r in rows})
code = {s: 's%02d' % i for i, s in enumerate(srcs)}
out = []
for r in rows:
 q = float(r['q']); b = float(r['b'])
 if not (b in (0.10, 0.30) or (b in (0.05, 0.20, 0.50) and q == 0.90)):
  continue
 rec = [code[r['src']], r['ds'], r['score'].replace('Synergy_',''), r['prior'],
        r['q'], r['b'], r['frac'], r['meth'], r['M'], r['MT'], r['n']]
 for k in KEEP:
  v = r.get(k, '')
  rec.append('' if v == '' else '%.4f' % float(v))
 out.append(','.join(rec))
hdr = ','.join(BASE + KEEP)
open('paper.csv','w').write(hdr + '\n' + '\n'.join(out) + '\n')
H = json.load(open('summary_head.json'))
H['srccode'] = {v: k for k, v in code.items()}
json.dump(H, open('paper_head.json','w'), separators=(',',':'))
print('slim rows', len(out), 'bytes', os.path.getsize('paper.csv'), 'head', os.path.getsize('paper_head.json'))
