"""POST HOC (L2): diagnose H4-HB 'formal' result; within-cohort chromosome-blocked CV; leave-one-MS-layer-out transfer matrix.
Not part of any frozen preregistration; frozen A4 results/verdicts are unchanged."""
import sys,os,json,numpy as np,pandas as pd
PH='./posthoc_L1L2'; sys.path.insert(0,PH)
os.environ['PH_TAG']='L2'; os.makedirs(PH+'/results_L2',exist_ok=True)
from c10common import *
import h4 as H4mod
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
SE='./seer_ext/'; A4='./a4/'
M=derive(load())
M=M.merge(pd.read_parquet(SE+'ms_ext_labels.parquet')[['a','ms_ext_primary']],on='a',how='left')
M=M.merge(pd.read_parquet(A4+'a4_cs_fields.parquet')[['a','ms_gs_primary','ms_hb_primary']],on='a',how='left').reset_index(drop=True)
Fall=H4mod.features(M)
LAY={'MSPA':'ms_primary','Seer':'ms_ext_primary','GS':'ms_gs_primary','HOLBAEK':'ms_hb_primary'}
MODS={'full':H4mod.FULL,'B0':H4mod.B0,'B1':H4mod.B1}
rng0=np.random.default_rng(SEED)
out={}
# ---------- (a) frozen A1 models ----------
tr=M[M.ms_ext_primary.notna()&(M.train==1)]; ytr=tr.ms_ext_primary.astype(int).values
frozen={}; coefrows=[]
for nm,cols in MODS.items():
    _,lr,(med,mu,sd)=H4mod.fit_pred(Fall.loc[tr.index],ytr,Fall.loc[tr.index[:2]],cols); frozen[nm]=(cols,lr,med,mu,sd)
    coefrows.append(dict(model=nm,term='(intercept)',coef_std=float(lr.intercept_[0])))
    for c_,b in zip(cols,lr.coef_[0]): coefrows.append(dict(model=nm,term=c_,coef_std=float(b),train_mean=float(mu[c_]),train_sd=float(sd[c_])))
out['frozen_coef']=coefrows; out['frozen_train_n']=int(len(tr)); out['frozen_train_pos']=int(ytr.sum())
def score(nm,idx):
    cols,lr,med,mu,sd=frozen[nm]; X=Fall.loc[idx,cols].fillna(med); return lr.predict_proba((X-mu)/sd)[:,1]
# univariate directions per layer
uni=[]
sets={'Seer A1 training set (odd chr+X)':tr.index,**{k:M.index[M[v].notna()] for k,v in LAY.items()}}
labcol={'Seer A1 training set (odd chr+X)':'ms_ext_primary',**LAY}
for sname,idx in sets.items():
    y=M.loc[idx,labcol[sname]].astype(int).values
    for c_ in H4mod.FULL:
        x=Fall.loc[idx,c_]; x=x.fillna(x.median()).astype(float)
        if x.nunique()<2: uni.append(dict(set=sname,term=c_,n=len(y),note='constant')); continue
        z=(x-x.mean())/x.std()
        try:
            f=sm.GLM(y,sm.add_constant(z.values),family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(M.loc[idx,'ga'])[0]})
            b,p=float(f.params[1]),float(f.pvalues[1])
        except Exception as e: b,p=np.nan,np.nan
        uni.append(dict(set=sname,term=c_,n=int(len(y)),n_pos=int(y.sum()),OR_per_SD=float(np.exp(b)),P_cluster=p,AUC_single=float(roc_auc_score(y,x))))
out['univariate']=uni
# B0 decomposition
b0=[]
for sname,idx in sets.items():
    T=M.loc[idx]; y=T[labcol[sname]].astype(int); s0=pd.Series(score('B0',idx),index=idx)
    for cat,g in T.groupby('conc_cat'):
        b0.append(dict(set=sname,conc_cat=cat,n=int(len(g)),supported=int(y.loc[g.index].sum()),rate=float(y.loc[g.index].mean()),frozen_B0_pred=float(s0.loc[g.index].iloc[0])))
    b0.append(dict(set=sname,conc_cat='AUC of frozen B0',n=int(len(T)),rate=float(roc_auc_score(y,s0))))
out['B0_decomposition']=b0
# ---------- helpers ----------
def boot_auc(y,preds,genes,B=1000,seed=SEED):
    rng=np.random.default_rng(seed); ug=np.unique(genes); gi={g:np.where(genes==g)[0] for g in ug}; res={k:[] for k in preds}; res['dfull_B0']=[]; res['dfull_B1']=[]
    for _ in range(B):
        ix=np.concatenate([gi[g] for g in rng.choice(ug,len(ug),replace=True)])
        if y[ix].min()==y[ix].max(): continue
        a={k:roc_auc_score(y[ix],p[ix]) for k,p in preds.items()}
        for k in preds: res[k].append(a[k])
        res['dfull_B0'].append(a['full']-a['B0']); res['dfull_B1'].append(a['full']-a['B1'])
    return {k:(float(np.percentile(v,2.5)),float(np.percentile(v,97.5))) for k,v in res.items()}
def summarize(y,preds,genes):
    o={f'auc_{k}':float(roc_auc_score(y,p)) for k,p in preds.items()}
    o['dAUC_full_B0']=o['auc_full']-o['auc_B0']; o['dAUC_full_B1']=o['auc_full']-o['auc_B1']
    o['ci']=boot_auc(y,preds,genes); o['delong_full_vs_B1']=H4mod.delong(y,preds['full'],preds['B1']); o['delong_full_vs_B0']=H4mod.delong(y,preds['full'],preds['B0'])
    return o
# ---------- (b) within-layer chromosome-blocked CV ----------
def cv_layer(lab,reps=20,k=5):
    T=M[M[lab].notna()]; idx=T.index.values; y=T[lab].astype(int).values; ch=T.chrom.values; genes=T.ga.values
    chroms=pd.Series(ch).value_counts()
    oof_sum={m:np.zeros(len(idx)) for m in MODS}; oof_cnt=np.zeros(len(idx)); rep_auc=[]
    rng=np.random.default_rng(SEED)
    for r in range(reps):
        order=list(chroms.index); rng.shuffle(order); order=sorted(order,key=lambda c_:-chroms[c_]+rng.random()*0.1)
        fold={}; sizes=np.zeros(k)
        for c_ in order: f=int(np.argmin(sizes+rng.random(k)*1e-3)); fold[c_]=f; sizes[f]+=chroms[c_]
        fid=np.array([fold[c_] for c_ in ch]); pr={m:np.full(len(idx),np.nan) for m in MODS}
        for f in range(k):
            te=fid==f; trn=~te
            if te.sum()==0 or len(np.unique(y[trn]))<2: continue
            for m,cols in MODS.items():
                pr[m][te]=H4mod.fit_pred(Fall.loc[idx[trn]],y[trn],Fall.loc[idx[te]],cols)[0]
        ok=~np.isnan(pr['full'])
        rep_auc.append({m:float(roc_auc_score(y[ok],pr[m][ok])) for m in MODS})
        for m in MODS: oof_sum[m][ok]+=pr[m][ok]
        oof_cnt[ok]+=1
    ok=oof_cnt>0; preds={m:oof_sum[m][ok]/oof_cnt[ok] for m in MODS}
    ra=pd.DataFrame(rep_auc)
    o=dict(n=int(len(idx)),n_pos=int(y.sum()),n_chrom=int(len(chroms)),reps=reps,k=k,
           rep_mean={m:float(ra[m].mean()) for m in MODS},rep_p2_5={m:float(ra[m].quantile(0.025)) for m in MODS},rep_p97_5={m:float(ra[m].quantile(0.975)) for m in MODS},
           rep_dAUC_full_B1_mean=float((ra.full-ra.B1).mean()),rep_dAUC_full_B0_mean=float((ra.full-ra.B0).mean()))
    o['oof_avg']=summarize(y[ok],preds,genes[ok])
    return o
out['cv']={nm:cv_layer(lab) for nm,lab in LAY.items()}
# ---------- (b) transfer matrix ----------
tm=[]
for A,la in LAY.items():
    TA=M[M[la].notna()]; yA=TA[la].astype(int).values
    fits={m:H4mod.fit_pred(Fall.loc[TA.index],yA,Fall.loc[TA.index[:2]],cols)[1:] for m,cols in MODS.items()}
    for Bn,lb in LAY.items():
        if A==Bn:
            cv=out['cv'][A]['oof_avg']; tm.append(dict(train=A,test=Bn,mode='within-layer CV (chromosome-blocked, 5-fold x20)',n_test=out['cv'][A]['n'],n_pos=out['cv'][A]['n_pos'],**{k:v for k,v in cv.items() if k.startswith(('auc','dAUC'))},ci=json.dumps(cv['ci']),P_full_vs_B1=cv['delong_full_vs_B1']['P_one_sided'])); continue
        TB=M[M[lb].notna()&~M.a.isin(TA.a)]; yB=TB[lb].astype(int).values
        if len(TB)<20 or yB.min()==yB.max(): tm.append(dict(train=A,test=Bn,mode='transfer',n_test=int(len(TB)),note='too few / one class')); continue
        preds={}
        for m,cols in MODS.items():
            lr,(med,mu,sd)=fits[m]; X=Fall.loc[TB.index,cols].fillna(med); preds[m]=lr.predict_proba((X-mu)/sd)[:,1]
        s=summarize(yB,preds,TB.ga.values)
        tm.append(dict(train=A,test=Bn,mode='transfer (test CS not in training set)',n_train=int(len(TA)),n_test=int(len(TB)),n_pos=int(yB.sum()),**{k:v for k,v in s.items() if k.startswith(('auc','dAUC'))},ci=json.dumps(s['ci']),P_full_vs_B1=s['delong_full_vs_B1']['P_one_sided']))
out['transfer']=tm
dump(out,'ph_L2.json'); print('done')
