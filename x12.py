import numpy as np,sys,json,os
exec(open('x1.py').read().split("DS=[('oneil")[0])
def alloc(rs,n,T):
 cnt=np.zeros(n,dtype=np.int64)
 if T>=n:
  cnt[:]=T//n;ex=int(T-cnt.sum())
  if ex>0: cnt[rs.choice(n,ex,replace=False)]+=1
 elif T>0: cnt[rs.choice(n,int(T),replace=False)]=1
 return cnt
def draw(rs,cnt,fv,sig):
 s=np.zeros(cnt.size);ok=cnt>0
 s[ok]=cnt[ok]*fv[ok]+sig*np.sqrt(cnt[ok])*rs.normal(0,1,int(ok.sum()))
 return s
def apt2(rs,fv,sig,T,tau,eps,nb=200):
 n=fv.size;cnt=np.zeros(n,dtype=np.int64)
 if T<=n:
  cnt=alloc(rs,n,T);return draw(rs,cnt,fv,sig),cnt
 s=fv+rs.normal(0,sig,n);cnt[:]=1;rem=int(T)-n;bs=max(1,rem//nb)
 while rem>0:
  b=int(min(bs,rem,n));pr=np.sqrt(cnt)*(np.abs(s/cnt-tau)+eps)
  sel=(np.arange(n) if b>=n else np.argpartition(pr,b-1)[:b])
  s[sel]+=fv[sel]+rs.normal(0,sig,b);cnt[sel]+=1;rem-=b
 return s,cnt
def post(prior,s,cnt,sig,sd):
 ok=cnt>0;yb=np.zeros(cnt.size);yb[ok]=s[ok]/cnt[ok]
 if ok.sum()>20:
  rr=float(np.mean((yb[ok]-prior[ok])**2)-sig**2*float(np.mean(1.0/cnt[ok])));vh=max(rr,1e-4*sd**2)
 else: vh=sd**2
 lam=sig**2/vh
 return np.where(ok,(lam*prior+cnt*yb)/(lam+cnt),prior)
DS=[('oneil','df_OncoPolyPharmacology.pkl','Y',8),('dcomb','df_DrugComb.pkl','Synergy_Loewe',12)]
FR=[0.0025,0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12]
BW=[0.05,0.10,0.20,0.30,0.50];QS=[0.80,0.85,0.90,0.95,0.99];SEEDS=5;A1=0.30
MET=['blind','passive','lse','band_inf']
which=sys.argv[1] if len(sys.argv)>1 else 'all'
for g,p,scol,r in DS:
 if which!='all' and which!=g: continue
 res=[]
 F=BLD(p,scol);v=np.isfinite(F);Fv=F[v];N=int(v.sum());P,m=F.shape;sd=float(Fv.std());sig=float(os.environ.get('SIGF','0.50'))*sd
 G=F.copy();G[~v]=Fv.mean();U,S,Vt=np.linalg.svd(G,full_matrices=False);Lv=((U[:,:r]*S[:r])@Vt[:r])[v];K=N
 TAU={q:float(np.quantile(Fv,q)) for q in QS}
 BAND={(q,b):np.where(np.abs(Lv-TAU[q])<=b*sd)[0] for q in QS for b in BW}
 print(json.dumps({'k':'meta','ds':g,'N':N,'P':P,'m':m,'r':r,'sd':sd,'sig':sig,
 'M':{'%.2f_%.2f'%(q,b):int(BAND[(q,b)].size) for q in QS for b in BW}}),flush=True)
 for frac in FR:
  T=max(int(frac*N),1)
  for seed in range(SEEDS):
   rs=np.random.default_rng(9176*int(frac*100000+1)+seed*7919)
   cn=np.bincount(rs.integers(0,K,T),minlength=K);sm=draw(rs,cn,Fv,sig)
   Ym=np.zeros((P,m));C=np.zeros((P,m));Ym[v]=np.where(cn>0,sm/np.maximum(cn,1),0.0);C[v]=cn
   Lp=fit(Ym,C,r)[v]
   T1=max(int(A1*T),1);c1=np.bincount(rs.integers(0,K,T1),minlength=K);s1=draw(rs,c1,Fv,sig)
   Y1=np.zeros((P,m));C1=np.zeros((P,m));Y1[v]=np.where(c1>0,s1/np.maximum(c1,1),0.0);C1[v]=c1
   Lh=fit(Y1,C1,r)[v]
   for q in QS:
    tau=TAU[q];Ht=Fv>=tau
    vb=post(np.full(K,float(Fv.mean())),sm,cn,sig,sd)>=tau
    vp=post(Lp,sm,cn,sig,sd)>=tau
    f1b=f1(vb,Ht);f1p=f1(vp,Ht)
    for b in BW:
     ar=BAND[(q,b)]
     if ar.size==0: continue
     def rec(mt,pred,ff):
      res.append({'ds':g,'q':q,'b':b,'frac':frac,'seed':seed,'meth':mt,'M':int(ar.size),
      'berr':float((pred[ar]!=Ht[ar]).mean()),'f1':ff})
     rec('blind',vb,f1b);rec('passive',vp,f1p)
     ci=np.zeros(K,dtype=np.int64);si=np.zeros(K);c2=alloc(rs,ar.size,T)
     ci[ar]=c2;si[ar]=draw(rs,c2,Fv[ar],sig);pi=post(Lv,si,ci,sig,sd)>=tau
     rec('band_inf',pi,f1(pi,Ht))
     arh=np.where(np.abs(Lh-tau)<=2*b*sd)[0];ct=c1.copy();st=s1.copy();rem=T-T1
     if arh.size>0 and rem>0:
      s2,c2b=apt2(rs,Fv[arh],sig,rem,tau,0.02*sd);ct[arh]+=c2b;st[arh]+=s2
     pl=post(Lh,st,ct,sig,sd)>=tau
     rec('lse',pl,f1(pl,Ht))
 json.dump(res,open('x12%s_%s.json'%(os.environ.get('TAG',''),g),'w'))
 for q in QS:
  for b in BW:
   for frac in FR:
    sub=[x for x in res if x['q']==q and x['b']==b and x['frac']==frac]
    if not sub: continue
    o={'k':'g','ds':g,'q':q,'b':b,'frac':frac,'M':sub[0]['M'],'NoverM':round(N/max(sub[0]['M'],1),3)}
    for mt in MET:
     e=[x['berr'] for x in sub if x['meth']==mt]
     o[mt]=round(float(np.mean(e)),5);o[mt+'_s']=round(float(np.std(e,ddof=1)),5)
    print(json.dumps(o),flush=True)
 print('DONE_'+g,flush=True)
print('DONE',flush=True)
