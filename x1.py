import numpy as np,pandas as pd,sys,json
def BLD(p,s):
    d=pd.read_pickle(p);a=d['Drug1_ID'].astype(str).values;b=d['Drug2_ID'].astype(str).values
    k=np.where(a<b,a+'|'+b,b+'|'+a)
    t=pd.DataFrame({'k':k,'c':d['Cell_Line_ID'].astype(str).values,'y':pd.to_numeric(d[s],errors='coerce').values})
    return t.pivot_table(index='k',columns='c',values='y',aggfunc='mean').values.astype(float)
def fit(Ym,C,r,lam=1.0,it=60):
    X=np.where(C>0,Ym,Ym[C>0].mean())
    for _ in range(it):
        U,S,Vt=np.linalg.svd(X,full_matrices=False);Xr=(U[:,:r]*S[:r])@Vt[:r]
        X=np.where(C>0,(C*Ym+lam*Xr)/(C+lam),Xr)
    U,S,Vt=np.linalg.svd(X,full_matrices=False);return (U[:,:r]*S[:r])@Vt[:r]
def f1(pred,true):
    tp=float((pred&true).sum());fp=float((pred&~true).sum());fn=float(((~pred)&true).sum())
    return 0.0 if tp==0 else 2*tp/(2*tp+fp+fn)
def complete(rs,K,T,i0,j0,P,m,Fv,sig,r):
    pk=rs.integers(0,K,T);Y=np.zeros((P,m));C=np.zeros((P,m))
    np.add.at(Y,(i0[pk],j0[pk]),Fv[pk]+rs.normal(0,sig,T));np.add.at(C,(i0[pk],j0[pk]),1.0)
    return fit(np.where(C>0,Y/np.maximum(C,1),0.0),C,r),Y,C
def unif_band(rs,ar,T,Fv,sig):
    n=ar.size;cnt=np.zeros(n,dtype=np.int64)
    if T>=n:
        cnt[:]=T//n;ex=T-int(cnt.sum())
        if ex>0: cnt[rs.choice(n,ex,replace=False)]+=1
    else:
        cnt[rs.choice(n,int(T),replace=False)]=1
    ok=cnt>0;mu=np.zeros(n)
    mu[ok]=Fv[ar[ok]]+sig/np.sqrt(cnt[ok])*rs.normal(0,1,int(ok.sum()))
    return mu,ok
def apt(rs,ar,Fv,sig,T,tau,eps,nb=200):
    n=ar.size;cnt=np.zeros(n);s=np.zeros(n)
    if T<=n:
        mu,ok=unif_band(rs,ar,T,Fv,sig);return mu,ok
    s=Fv[ar]+rs.normal(0,sig,n);cnt[:]=1.0;rem=int(T)-n;bs=max(1,rem//nb)
    while rem>0:
        b=int(min(bs,rem,n));pr=np.sqrt(cnt)*(np.abs(s/cnt-tau)+eps)
        sel=(np.arange(n) if b>=n else np.argpartition(pr,b-1)[:b])
        s[sel]+=Fv[ar[sel]]+rs.normal(0,sig,b);cnt[sel]+=1.0;rem-=b
    return s/cnt,np.ones(n,dtype=bool)
def lev(L,r,P,m):
    U,S,Vt=np.linalg.svd(L,full_matrices=False);u=U[:,:r];w=Vt[:r,:].T
    return P*(u**2).sum(1)/r,m*(w**2).sum(1)/r
DS=[('oneil','df_OncoPolyPharmacology.pkl','Y',8),('dcomb','df_DrugComb.pkl','Synergy_Loewe',12)]
QS=[0.80,0.85,0.90,0.95,0.99];FR=[0.05,0.10,0.20,0.30,0.50,0.80];SEEDS=5;ASPLIT=0.30
METH=['blind','passive','band_inf','struct_na','lse']
which=sys.argv[1] if len(sys.argv)>1 else 'all'
res=[]
for g,p,scol,r in DS:
    if which!='all' and which!=g: continue
    F=BLD(p,scol);v=np.isfinite(F);Fv=F[v];N=int(v.sum());P,m=F.shape;sd=float(Fv.std())
    sig=0.50*sd;bw=0.30*sd
    G=F.copy();G[~v]=Fv.mean();U,S,Vt=np.linalg.svd(G,full_matrices=False);Lm=(U[:,:r]*S[:r])@Vt[:r]
    Lv=Lm[v];ix=np.argwhere(v);i0=ix[:,0];j0=ix[:,1];K=N
    mu_i,nu_j=lev(Lm,r,P,m)
    eta=Fv-Lv
    print(json.dumps({'k':'meta','ds':g,'N':N,'P':P,'m':m,'r':r,'sd':sd,'sig':sig,'bw':bw,
      'mu_max':float(mu_i.max()),'nu_max':float(nu_j.max()),'mu_med':float(np.median(mu_i)),
      'eta_max_sd':float(np.abs(eta).max()/sd),'eta_l2_over_N':float((eta**2).mean()),
      'ev':float(1-(eta**2).sum()/((Fv-Fv.mean())**2).sum()),'mu2r':float(mu_i.max()**2*r*(P+m)),'rPm':float(r*(P+m))}),flush=True)
    for q in QS:
        tau=float(np.quantile(Fv,q));Ht=Fv>=tau
        Bt=np.abs(Lv-tau)<=bw;M=int(Bt.sum());art=np.where(Bt)[0]
        rows_b=np.unique(i0[art]);cols_b=np.unique(j0[art])
        muB=float(mu_i[rows_b].max()) if rows_b.size else 0.0
        nuB=float(nu_j[cols_b].max()) if cols_b.size else 0.0
        pB=float((np.abs(eta)>bw).mean());l2B=float((eta**2).mean()/bw**2)
        print(json.dumps({'k':'inst','ds':g,'q':q,'tau':tau,'M':M,'NoverM':N/max(M,1),'pB':pB,
          'l2_over_B2':l2B,'muB':muB,'nuB':nuB,'muB2r':float(muB**2*r*(P+m))}),flush=True)
        for frac in FR:
            T=int(frac*N)
            for meth in METH:
                for seed in range(SEEDS):
                    rs=np.random.default_rng(abs(hash((g,q,frac,meth,seed)))%(2**31) if False else (int(q*100)*1000003+int(frac*1000)*10007+seed*101+len(g)*7+METH.index(meth)*13))
                    if meth=='blind':
                        pk=rs.integers(0,K,T);cnt=np.bincount(pk,minlength=K).astype(float)
                        ok=cnt>0;mu=np.zeros(K)
                        mu[ok]=Fv[ok]+sig/np.sqrt(cnt[ok])*rs.normal(0,1,int(ok.sum()))
                        pred=np.zeros(K,dtype=bool);pred[ok]=mu[ok]>=tau;Me=int(ok.sum())
                    elif meth=='passive':
                        Lh,_,_=complete(rs,K,T,i0,j0,P,m,Fv,sig,r);pred=Lh[v]>=tau;Me=0
                    elif meth=='band_inf':
                        pred=Lv>=tau;Me=M
                        if art.size>0:
                            mu,ok=unif_band(rs,art,T,Fv,sig);pr=pred[art].copy();pr[ok]=mu[ok]>=tau;pred[art]=pr
                    else:
                        T1=max(int(ASPLIT*T),1);Lh,_,_=complete(rs,K,T1,i0,j0,P,m,Fv,sig,r);Lhv=Lh[v]
                        arh=np.where(np.abs(Lhv-tau)<=2*bw)[0];pred=Lhv>=tau;Me=int(arh.size);rem=T-T1
                        if arh.size>0 and rem>0:
                            if meth=='struct_na': mu,ok=unif_band(rs,arh,rem,Fv,sig)
                            else: mu,ok=apt(rs,arh,Fv,sig,rem,tau,0.02*sd)
                            pr=pred[arh].copy();pr[ok]=mu[ok]>=tau;pred[arh]=pr
                    res.append({'k':'run','ds':g,'q':q,'frac':frac,'meth':meth,'seed':seed,'T':T,
                      'f1':f1(pred,Ht),'err':float((pred!=Ht).mean()),'Mband':Me})
            print(json.dumps({'k':'agg','ds':g,'q':q,'frac':frac,'M':M,'NoverM':N/max(M,1),
              **{mt+'_f1':float(np.mean([x['f1'] for x in res if x['meth']==mt and x['q']==q and x['frac']==frac and x['ds']==g])) for mt in METH},
              **{mt+'_f1sd':float(np.std([x['f1'] for x in res if x['meth']==mt and x['q']==q and x['frac']==frac and x['ds']==g])) for mt in METH},
              **{mt+'_err':float(np.mean([x['err'] for x in res if x['meth']==mt and x['q']==q and x['frac']==frac and x['ds']==g])) for mt in METH}}),flush=True)
json.dump(res,open('x1_%s.json'%which,'w'))
print('DONE',flush=True)
