"""Emit every figure, table body and inline number the paper needs."""
import csv, json, os, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.makedirs('Figures', exist_ok=True)
NUMF = {'q','b','frac','n','berrT','berrT_lo','berrT_hi','berr','berrBal','f1','f1_lo','f1_hi',
        'posfrac','bfp','bfn','cost3','cost10','gcost3','MT','Mhat','pri_e','prip_e','prih_e',
        'bprec','brec','nsel','lam','gerr','prec','rec','gfp','gfn','gcost10','M','pred'}
R = []
for d in csv.DictReader(open('summary_y.csv')):
    r = dict(d)
    for k in NUMF:
        v = d.get(k, '')
        r[k] = None if v in ('', None) else float(v)
    R.append(r)
IDX = {}
for r in R:
    IDX[(r['srcname'], r['meth'], round(r['q'],2), round(r['b'],2), round(r['frac'],6))] = r

SRC = {('oneil','svd','homo'):'xy_oneil_Y_svd_homo',
       ('dcomb','svd','homo'):'xy_dcomb_Synergy_Loewe_svd_homo',
       ('oneil','svd','emp'):'xy_oneil_Y_svd_emp',
       ('dcomb','svd','emp'):'xy_dcomb_Synergy_Loewe_svd_emp',
       ('oneil','svd','ig'):'xy_oneil_Y_svd_homo_ig',
       ('dcomb','svd','ig'):'xy_dcomb_Synergy_Loewe_svd_homo_ig',
       ('oneil','svd','sp'):'xy_oneil_Y_svd_homo_sp',
       ('dcomb','svd','sp'):'xy_dcomb_Synergy_Loewe_svd_homo_sp',
       ('oneil','gbt','homo'):'xy_oneil_Y_gbt_homo',
       ('dcomb','gbt','homo'):'xy_dcomb_Synergy_Loewe_gbt_homo',
       ('oneil','clf','homo'):'xy_oneil_Y_clf_homo',
       ('dcomb','clf','homo'):'xy_dcomb_Synergy_Loewe_clf_homo',
       ('oneil','clfu','homo'):'xy_oneil_Y_clfu_homo',
       ('oneil','r2','ab'):'xy_oneil_Y_svd_homo_r2',
       ('oneil','r5','ab'):'xy_oneil_Y_svd_homo_r5',
       ('oneil','r20','ab'):'xy_oneil_Y_svd_homo_r20',
       ('oneil','r32','ab'):'xy_oneil_Y_svd_homo_r32',
       ('oneil','a0.10','ab'):'xy_oneil_Y_svd_homo_a0.10',
       ('oneil','a0.50','ab'):'xy_oneil_Y_svd_homo_a0.50',
       ('dcomb','r4','ab'):'xy_dcomb_Synergy_Loewe_svd_homo_r4',
       ('dcomb','r20','ab'):'xy_dcomb_Synergy_Loewe_svd_homo_r20',
       ('dcomb','r32','ab'):'xy_dcomb_Synergy_Loewe_svd_homo_r32',
       ('dcomb','clfu','homo'):'xy_dcomb_Synergy_Loewe_clfu_homo'}
NN = {'oneil':22737,'dcomb':287268}
MM = {'oneil':39,'dcomb':59}
RR = {'oneil':8,'dcomb':12}
DSN = {'oneil':"O'Neil",'dcomb':'DrugComb'}
Q0, B0 = 0.90, 0.30

def g(key, meth, metric, frac, q=Q0, b=B0):
    # INFOGAIN was run into its own file; route it there transparently
    if meth == 'ig' and key[2] in ('homo',):
        key = (key[0], key[1], 'ig')
    if key not in SRC:
        return None
    r = IDX.get((SRC[key], meth, round(q,2), round(b,2), round(frac,6)))
    return None if r is None else r.get(metric)

def has(key):
    return SRC.get(key) and any(k[0]==SRC[key] for k in IDX)

def f3(x): return '--' if x is None else format(round(x+1e-12,3),'.3f')
def com(n): return format(int(round(n)),',').replace(',','{,}')

FR = sorted({r['frac'] for r in R if r['srcname']==SRC[('oneil','svd','homo')]})
BW = sorted({r['b'] for r in R if r['srcname']==SRC[('oneil','svd','homo')]})
QS = sorted({r['q'] for r in R if r['srcname']==SRC[('oneil','svd','homo')]})
T3 = [x for x in (0.04,0.32,1.28) if x in FR]

LBL = {'null0':'NONE','null1':'ALL','meas':'MEAS','blind':'UNIF','passive':'STRUCT',
       'band_hat':'BAND','lse':'ADAPT','ig':'INFOGAIN','band_inf':'ORACLE',
       'blind_split':'UNIF-S','passive_split':'STRUCT-S'}
ORD = ['null0','null1','meas','blind','passive','band_hat','lse','ig','band_inf']
TARGET = ['band_hat','lse','ig','band_inf']

facts = {'grid': {'fracs': FR, 'bands': BW, 'quants': QS}}

# ---------------------------------------------------------------- leak
leak = {}
for ds in ('oneil','dcomb'):
    key = (ds,'svd','homo')
    if not has(key): continue
    leak[ds] = {}
    for b in BW:
        row = {}
        ctrl = None
        for m in ORD + ['blind']:
            e_s = g(key, m, 'berr', 0.32, b=b)      # model-defined band
            e_t = g(key, m, 'berrT', 0.32, b=b)     # model-free band
            if e_s is None or e_t is None or e_t == 0: continue
            row[m] = {'B': round(e_s,5), 'Bstar': round(e_t,5), 'ratio': round(e_s/e_t,4)}
        if 'blind' in row: ctrl = row['blind']['ratio']
        for m in row:
            row[m]['inflation'] = round(ctrl/row[m]['ratio'],4) if ctrl and row[m]['ratio'] else None
        leak[ds][b] = row
facts['leak'] = leak

# ---------------------------------------------------------------- headline facts
hl = {}
for ds in ('oneil','dcomb'):
    key = (ds,'svd','homo')
    if not has(key): continue
    hl[ds] = {}
    for fr in T3:
        hl[ds][fr] = {m: {k: g(key,m,k,fr) for k in
                          ('berrT','berrBal','berr','f1','bfp','bfn','cost3','cost10','posfrac','MT')}
                      for m in ORD}
facts['headline'] = hl
facts['bands'] = {ds: {b: {'Mstar': g((ds,'svd','homo'),'blind','MT',0.32,b=b),
                           'posfrac': g((ds,'svd','homo'),'blind','posfrac',0.32,b=b)}
                       for b in BW} for ds in ('oneil','dcomb') if has((ds,'svd','homo'))}
for ds, dd in facts['bands'].items():
    for b, o in dd.items():
        if o['Mstar']:
            o['NoverM'] = round(NN[ds]/o['Mstar'],2)
            o['deff'] = round(o['Mstar']*(RR[ds]-1)/MM[ds],1)

# ---------------------------------------------------------------- priors
pri = {}
for ds in ('oneil','dcomb'):
    pri[ds] = {}
    for pk, nm in (('svd','spine'),('gbt','boosting'),('clf','sign'),('clfu','sign-g')):
        key = (ds,pk,'homo')
        if not has(key): continue
        pri[ds][nm] = {fr: {m: {'berrT': g(key,m,'berrT',fr), 'f1': g(key,m,'f1',fr)}
                            for m in ('null0','blind','passive','band_hat','lse','band_inf')}
                       for fr in T3}
facts['priors'] = pri

# ---------------------------------------------------------------- noise
noise = {}
for ds in ('oneil','dcomb'):
    noise[ds] = {}
    for nk in ('homo','emp'):
        key = (ds,'svd',nk)
        if not has(key): continue
        noise[ds][nk] = {fr: {m: g(key,m,'berrT',fr) for m in ORD} for fr in T3}
facts['noise'] = noise

# ---------------------------------------------------------------- split
split = {}
for ds in ('oneil','dcomb'):
    key = (ds,'svd','sp')
    if not has(key): continue
    split[ds] = {fr: {m: {'berrT': g(key,m,'berrT',fr), 'pred': g(key,m,'pred',fr)}
                      for m in ('null0','blind','passive','blind_split','passive_split')}
                 for fr in T3}
facts['split'] = split

abl = {}
ABK = {'oneil': [('r2','rank 2'),('r5','rank 5'),('r20','rank 20'),('r32','rank 32'),
                 ('a0.10','split 0.10'),('a0.50','split 0.50')],
       'dcomb': [('r4','rank 4'),('r20','rank 20'),('r32','rank 32')]}
for _ds, _lst in ABK.items():
    abl[_ds] = {}
    for _tg, _lab in _lst:
        _k = (_ds, _tg, 'ab')
        if not has(_k):
            continue
        abl[_ds][_lab] = {fr: {m: g(_k, m, 'berrT', fr)
                               for m in ('null0','blind','passive','band_hat','lse')}
                          for fr in T3}
facts['ablate'] = abl
print(json.dumps({'k':'ablate','have':{k:list(v) for k,v in abl.items()}}), flush=True)

json.dump(facts, open('facts_y.json','w'), indent=1, default=str)

# ---------------------------------------------------------------- figure data
fig = {'fracs': FR, 'bands': BW, 'quants': QS, 'series': {}}
for ds in ('oneil', 'dcomb'):
    for nk in ('homo', 'emp', 'ig', 'sp'):
        key = (ds, 'svd', nk)
        if not has(key):
            continue
        for m in ORD + ['blind_split', 'passive_split']:
            d = {}
            for metric in ('berrT', 'berrT_lo', 'berrT_hi', 'f1', 'f1_lo', 'f1_hi',
                           'cost3', 'berrBal', 'pred'):
                vals = [g(key, m, metric, fr) for fr in FR]
                if any(v is not None for v in vals):
                    d[metric] = vals
            if d:
                fig['series']['%s|%s|%s' % (ds, nk, m)] = d
    for pk in ('gbt', 'clf', 'clfu'):
        key = (ds, pk, 'homo')
        if not has(key):
            continue
        for m in ORD:
            d = {}
            for metric in ('berrT', 'f1', 'cost3'):
                vals = [g(key, m, metric, fr) for fr in FR]
                if any(v is not None for v in vals):
                    d[metric] = vals
            if d:
                fig['series']['%s|%s|%s' % (ds, pk, m)] = d
fig['leak'] = leak
fig['bands_info'] = facts['bands']
json.dump(fig, open('figdata.json', 'w'), separators=(',', ':'), default=str)
print(json.dumps({'k': 'figdata', 'series': len(fig['series'])}), flush=True)
print(json.dumps({'k':'facts','leak_ds':list(leak),'priors':{k:list(v) for k,v in pri.items()},
                  'noise':{k:list(v) for k,v in noise.items()},'split':list(split)}), flush=True)

# ================================================================ FIGURES
plt.rcParams.update({
    'font.family':'serif','font.serif':['Times','DejaVu Serif'],'font.size':8,
    'axes.labelsize':8,'axes.titlesize':8,'legend.fontsize':7,'xtick.labelsize':7,
    'ytick.labelsize':7,'axes.linewidth':0.6,'grid.linewidth':0.4,'lines.linewidth':1.3,
    'savefig.bbox':'tight','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
COL = {'null0':'#999999','null1':'#cccccc','meas':'#56B4E9','blind':'#000000',
       'passive':'#D55E00','band_hat':'#0072B2','lse':'#009E73','ig':'#CC79A7',
       'band_inf':'#E69F00'}
STY = {'null0':(0,(1,2)),'null1':(0,(1,4)),'meas':(0,(3,1,1,1)),'blind':(0,(1,1)),
       'passive':'-','band_hat':'-','lse':'-','ig':'-','band_inf':(0,(4,2))}
MRK = {'null0':None,'null1':None,'meas':'P','blind':'o','passive':'s','band_hat':'^',
       'lse':'D','ig':'*','band_inf':'v'}
def grid(ax):
    ax.grid(True, color='#cccccc', alpha=0.6); ax.set_axisbelow(True)

# ---- Figure: leakage inflation vs band width
def fig_leak():
    fg, axs = plt.subplots(1, 2, figsize=(6.4, 2.3), sharey=True)
    for ax, ds in zip(axs, ('oneil','dcomb')):
        if ds not in leak: continue
        for m in TARGET:
            xs, ys = [], []
            for b in BW:
                o = leak[ds].get(b, {}).get(m)
                if o and o.get('inflation'): xs.append(b); ys.append(o['inflation'])
            if len(xs) > 1:
                ax.plot(xs, ys, marker=MRK[m], markersize=3.5, markeredgewidth=0,
                        color=COL[m], label=LBL[m])
        ax.axhline(1.0, color='#000000', linewidth=0.9)
        ax.text(BW[-1], 1.03, 'non-targeting control', ha='right', va='bottom', fontsize=6.5)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xticks(BW); ax.set_xticklabels([str(z) for z in BW])
        ax.set_xlabel(r'band half-width $b$ ($\times\sigma_F$)')
        ax.set_title(DSN[ds], loc='left')
        grid(ax)
    axs[0].set_ylabel('apparent gain inflation')
    axs[0].legend(frameon=False, loc='upper right')
    fg.tight_layout(); fg.savefig('Figures/leak.pdf'); plt.close(fg)
    print('  Figures/leak.pdf', flush=True)

# ---- Figure: headline, F1 and boundary error vs budget
def fig_headline():
    fg, axs = plt.subplots(2, 2, figsize=(6.4, 3.6), sharex=True)
    hh, ll = [], []
    for col, ds in enumerate(('oneil','dcomb')):
        key = (ds,'svd','homo')
        if not has(key): continue
        for row, (metric, ylab) in enumerate(((('f1'),'global $F_1$'), ('berrT','boundary error'))):
            ax = axs[row][col]
            for m in ORD:
                if m == 'null1' and row == 0: continue
                xs = [fr for fr in FR if g(key,m,metric,fr) is not None]
                if not xs: continue
                ys = [g(key,m,metric,fr) for fr in xs]
                if m in ('null0','null1'):
                    ax.axhline(ys[0], color=COL[m], linestyle=STY[m], linewidth=1.1)
                    if row == 1 and col == 0:
                        hh.append(plt.Line2D([],[],color=COL[m],linestyle=STY[m])); ll.append(LBL[m])
                    continue
                lo = [g(key,m,metric+'_lo',fr) or ys[i] for i,fr in enumerate(xs)]
                hi = [g(key,m,metric+'_hi',fr) or ys[i] for i,fr in enumerate(xs)]
                ax.fill_between(xs, lo, hi, color=COL[m], alpha=0.13, linewidth=0)
                ln, = ax.plot(xs, ys, linestyle=STY[m], color=COL[m], marker=MRK[m],
                              markersize=3.0, markeredgewidth=0)
                if row == 0 and col == 0: hh.append(ln); ll.append(LBL[m])
            ax.set_xscale('log')
            if row == 1: ax.set_xlabel(r'query budget $T/N$')
            if col == 0: ax.set_ylabel(ylab)
            if row == 0: ax.set_title(DSN[ds], loc='left')
            grid(ax)
    fg.tight_layout(rect=(0,0.10,1,1))
    fg.legend(hh, ll, frameon=False, loc='lower center', ncol=5, handlelength=2.0,
              columnspacing=1.3, bbox_to_anchor=(0.5,-0.02))
    fg.savefig('Figures/headline.pdf'); plt.close(fg)
    print('  Figures/headline.pdf', flush=True)

# ---- Figure: prior class trained for the sign
def fig_sign():
    names = [('spine','svd'),('boosting','gbt'),('sign','clf'),('sign-g','clfu')]
    fg, axs = plt.subplots(1, 2, figsize=(6.4, 2.3))
    w = 0.19
    for ax, ds in zip(axs, ('oneil','dcomb')):
        avail = [(nm,pk) for nm,pk in names if has((ds,pk,'homo'))]
        if not avail: continue
        x = np.arange(len(T3))
        for i,(nm,pk) in enumerate(avail):
            ys = [g((ds,pk,'homo'),'passive','berrT',fr) for fr in T3]
            ys = [np.nan if y is None else y for y in ys]
            ax.bar(x + (i-len(avail)/2+0.5)*w, ys, w, label=nm,
                   color=['#D55E00','#0072B2','#009E73','#CC79A7'][i], edgecolor='none')
        u = [g((ds,'svd','homo'),'blind','berrT',fr) for fr in T3]
        o = [g((ds,'svd','homo'),'band_inf','berrT',fr) for fr in T3]
        n0 = [g((ds,'svd','homo'),'null0','berrT',fr) for fr in T3]
        for arr,st,lb,c in ((u,(0,(1,1)),'UNIF','#000000'),(o,(0,(4,2)),'ORACLE','#E69F00'),
                            (n0,(0,(1,2)),'NONE','#999999')):
            if all(z is not None for z in arr):
                ax.plot(x, arr, linestyle=st, color=c, marker='', linewidth=1.1, label=lb)
        ax.set_xticks(x); ax.set_xticklabels([('%gN'%fr) for fr in T3])
        ax.set_xlabel(r'query budget $T$'); ax.set_title(DSN[ds], loc='left')
        grid(ax)
    axs[0].set_ylabel('boundary error')
    axs[0].legend(frameon=False, loc='best', fontsize=6, ncol=2)
    fg.tight_layout(); fg.savefig('Figures/sign.pdf'); plt.close(fg)
    print('  Figures/sign.pdf', flush=True)

# ---- Figure: calibrated vs optimistic noise
def fig_noise():
    fg, axs = plt.subplots(1, 2, figsize=(6.4, 2.3), sharey=False)
    for ax, ds in zip(axs, ('oneil','dcomb')):
        for nk, st, tag in (('homo','-',r'$\sigma=0.50\sigma_F$, homoskedastic'),
                            ('emp',(0,(4,2)),'calibrated, heteroskedastic')):
            key = (ds,'svd',nk)
            if not has(key): continue
            for m in ('blind','passive','lse','null0'):
                xs = [fr for fr in FR if g(key,m,'berrT',fr) is not None]
                if not xs: continue
                ys = [g(key,m,'berrT',fr) for fr in xs]
                ax.plot(xs, ys, linestyle=st, color=COL[m], marker=MRK[m] if nk=='homo' else None,
                        markersize=2.6, markeredgewidth=0,
                        label=('%s, %s' % (LBL[m], 'optimistic' if nk=='homo' else 'calibrated')))
        ax.set_xscale('log'); ax.set_xlabel(r'query budget $T/N$')
        ax.set_title(DSN[ds], loc='left'); grid(ax)
    axs[0].set_ylabel('boundary error')
    axs[0].legend(frameon=False, loc='best', fontsize=5.5, ncol=2)
    fg.tight_layout(); fg.savefig('Figures/noise.pdf'); plt.close(fg)
    print('  Figures/noise.pdf', flush=True)

# ---- Figure: mechanism (identity)
def fig_mechanism():
    from scipy.stats import norm
    sig, lam = 1.0, 1.0
    x = np.linspace(-2.0, 2.0, 601); D = 0.35
    fg, ax = plt.subplots(figsize=(3.1, 2.15))
    cols = ['#0072B2','#009E73','#E69F00','#D55E00']
    for i, c in enumerate([1,2,4,8]):
        w = lam/(lam+c)
        ax.plot(x, norm.cdf(-(lam+c)/(sig*np.sqrt(c))*(D + w*(x-D))), color=cols[i], label='$c=%d$'%c)
    ax.axvline(0.0, color='#000000', linewidth=0.9)
    ax.axhline(0.5, color='#888888', linewidth=0.6, linestyle=(0,(2,2)))
    ax.axvspan(-2.0, 0.0, color='#D55E00', alpha=0.07, linewidth=0)
    ax.text(-1.02, 0.60, 'prior on the\nwrong side', ha='center', va='center', fontsize=6.5, color='#8a3800')
    ax.text(1.30, 0.62, 'prior on the\ncorrect side', ha='center', va='center', fontsize=6.5, color='#00506b')
    ax.set_xlabel(r'signed prior position $\mathrm{sign}(\Delta_e)(L_e-\tau)/\sigma_F$')
    ax.set_ylabel('error probability'); ax.set_ylim(0,1.0)
    ax.legend(frameon=False, loc='upper right', ncol=2, columnspacing=1.0)
    grid(ax); fg.tight_layout(); fg.savefig('Figures/mechanism.pdf'); plt.close(fg)
    print('  Figures/mechanism.pdf', flush=True)

# ---- Figure: asymmetric cost
def fig_cost():
    fg, axs = plt.subplots(1, 2, figsize=(6.4, 2.3))
    for ax, ds in zip(axs, ('oneil','dcomb')):
        key = (ds,'svd','homo')
        if not has(key): continue
        for m in ('null0','blind','passive','lse','ig','band_inf'):
            xs = [fr for fr in FR if g(key,m,'cost3',fr) is not None]
            if not xs: continue
            ys = [g(key,m,'cost3',fr) for fr in xs]
            if m == 'null0':
                ax.axhline(ys[0], color=COL[m], linestyle=STY[m], linewidth=1.1, label=LBL[m]); continue
            ax.plot(xs, ys, linestyle=STY[m], color=COL[m], marker=MRK[m], markersize=3.0,
                    markeredgewidth=0, label=LBL[m])
        ax.set_xscale('log'); ax.set_xlabel(r'query budget $T/N$')
        ax.set_title(DSN[ds], loc='left'); grid(ax)
    axs[0].set_ylabel(r'boundary cost $C_3$')
    axs[0].legend(frameon=False, loc='best', fontsize=6, ncol=2)
    fg.tight_layout(); fg.savefig('Figures/cost.pdf'); plt.close(fg)
    print('  Figures/cost.pdf', flush=True)

for f in (fig_leak, fig_headline, fig_sign, fig_noise, fig_mechanism, fig_cost):
    try: f()
    except Exception as e: print('  FIG FAIL', f.__name__, repr(e)[:160], flush=True)

# ================================================================ TABLES
out = []
NLJ = chr(10)

def emit(name, body):
    out.append('%% ---- %s\n%s' % (name, body))

# leak table
for ds in ('oneil','dcomb'):
    if ds not in leak: continue
    L = []
    for m in ORD:
        o = leak[ds].get(B0, {}).get(m)
        o10 = leak[ds].get(0.10, {}).get(m)
        if not o: continue
        L.append('\\textsc{%s} & %s & %s & %s & %s & %s & %s\\\\' %
                 (LBL[m].lower(), f3(o['Bstar']), f3(o['B']), f3(o['ratio']),
                  f3(o10['Bstar']) if o10 else '--', f3(o10['B']) if o10 else '--',
                  f3(o10['ratio']) if o10 else '--'))
    emit('tab:leak %s' % ds, '\n'.join(L))

# headline table
for metric in ('berrT','f1'):
    L = []
    for m in ORD:
        cs = []
        for ds in ('oneil','dcomb'):
            for fr in T3:
                cs.append(f3(g((ds,'svd','homo'), m, metric, fr)))
        L.append('\\textsc{%s} & %s\\\\' % (LBL[m].lower(), ' & '.join(cs)))
    emit('tab:headline %s' % metric, '\n'.join(L))

# cost table
L = []
for m in ORD:
    cs = []
    for ds in ('oneil','dcomb'):
        for k in ('bfp','bfn','cost3','cost10'):
            cs.append(f3(g((ds,'svd','homo'), m, k, 0.32)))
    L.append('\\textsc{%s} & %s\\\\' % (LBL[m].lower(), ' & '.join(cs)))
emit('tab:cost', '\n'.join(L))

# bands table
L = []
for ds in ('oneil','dcomb'):
    cs = []
    for b in BW:
        o = facts['bands'].get(ds, {}).get(b)
        if not o or not o.get('Mstar'): cs.append('--'); continue
        cs.append('%s / %s' % (com(o['Mstar']), com(o['deff'])))
    L.append('%s & %s\\\\' % (DSN[ds], ' & '.join(cs)))
emit('tab:bands', '\n'.join(L))

# prior table
L = []
for ds in ('oneil','dcomb'):
    for nm in ('spine','boosting','sign','sign-g'):
        if nm not in pri.get(ds, {}): continue
        cs = []
        for fr in T3:
            v = pri[ds][nm][fr]['passive']
            cs += [f3(v['berrT']), f3(v['f1'])]
        L.append('%s & %s & %s\\\\' % (DSN[ds], nm, ' & '.join(cs)))
emit('tab:priors', '\n'.join(L))

# noise table
L = []
for ds in ('oneil','dcomb'):
    for nk, lab in (('homo','optimistic'),('emp','calibrated')):
        if nk not in noise.get(ds, {}): continue
        cs = []
        for fr in T3:
            for m in ('null0','blind','passive','lse'):
                cs.append(f3(noise[ds][nk][fr].get(m)))
        L.append('%s & %s & %s\\\\' % (DSN[ds], lab, ' & '.join(cs)))
emit('tab:noise', '\n'.join(L))

# split table
L = []
for ds in ('oneil','dcomb'):
    if ds not in split: continue
    for m in ('blind','passive','blind_split','passive_split'):
        cs = []
        for fr in T3:
            v = split[ds][fr].get(m, {})
            cs += [f3(v.get('berrT')), f3(v.get('pred'))]
        L.append('%s & \\textsc{%s} & %s\\\\' % (DSN[ds], LBL[m].lower(), ' & '.join(cs)))
emit('tab:split', '\n'.join(L))

# q sweep
L = []
for ds in ('oneil','dcomb'):
    for m in ('null0','blind','passive','lse'):
        cs = [f3(g((ds,'svd','homo'), m, 'berrT', 0.32, q=q)) for q in QS]
        L.append('%s & \\textsc{%s} & %s\\\\' % (DSN[ds], LBL[m].lower(), ' & '.join(cs)))
emit('tab:qsweep', '\n'.join(L))

L = []
for ds in ('oneil','dcomb'):
    for lab, dd in facts.get('ablate', {}).get(ds, {}).items():
        cs = []
        for fr in T3:
            cs += [f3(dd[fr].get('blind')), f3(dd[fr].get('passive')), f3(dd[fr].get('lse'))]
        L.append('%s & %s & %s\\\\' % (DSN[ds], lab, ' & '.join(cs)))
emit('tab:ablate', chr(10).join(L))

L = []
for m in ORD:
    cs = []
    for ds in ('oneil','dcomb'):
        for fr in T3:
            mu = g((ds,'svd','homo'), m, 'berrT', fr)
            lo = g((ds,'svd','homo'), m, 'berrT_lo', fr)
            hi = g((ds,'svd','homo'), m, 'berrT_hi', fr)
            if mu is None:
                cs.append('--')
            elif lo is None or hi is None:
                cs.append(f3(mu))
            else:
                cs.append(r'%s$\\,\pm\\,$%s' % (f3(mu), f3((hi-lo)/2.0)))
    L.append('\\textsc{%s} & %s\\\\' % (LBL[m].lower(), ' & '.join(cs)))
emit('tab:ci berrT', NLJ.join(L))

L = []
for m in ORD:
    cs = []
    for ds in ('oneil','dcomb'):
        for fr in T3:
            cs.append(f3(g((ds,'svd','homo'), m, 'berrBal', fr)))
    L.append('\\textsc{%s} & %s\\\\' % (LBL[m].lower(), ' & '.join(cs)))
emit('tab:balanced', NLJ.join(L))

open('tables_y.txt','w').write('\n\n'.join(out) + '\n')
print(json.dumps({'k':'tables','blocks':len(out),'file':'tables_y.txt'}), flush=True)
print('MKO_DONE', flush=True)
