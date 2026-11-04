"""POST HOC (L1c): A1 + A4 MS endpoints with/without chrX, using frozen labels (no relabelling) and the frozen A1 model."""
import sys,os,json,numpy as np,pandas as pd
PH='./posthoc_L1L2'; sys.path.insert(0,PH)
from c10common import *
import h3 as H3mod, h4 as H4mod
SE='./seer_ext/'; A4='./a4/'
M=derive(load())
M=M.merge(pd.read_parquet(SE+'ms_ext_labels.parquet')[['a','ms_ext_primary','seer_ext_studies']],on='a',how='left')
F=pd.read_parquet(A4+'a4_cs_fields.parquet')[['a','ms_gs_primary','ms_hb_primary','ms_independent_combined']]; M=M.merge(F,on='a',how='left')
res={'n_cs':int(len(M)),'n_chrX':int((M.chrom=='X').sum())}
for lab in ['ms_primary','ms_ext_primary','ms_gs_primary','ms_hb_primary','ms_independent_combined']:
    res[f'n_det_{lab}']=int(M[lab].notna().sum()); res[f'n_det_chrX_{lab}']=int((M[lab].notna()&(M.chrom=='X')).sum())
res['H3_mspa'],_=H3mod.h3_block(M,'ms_primary')
res['H3_ext'],_=H3mod.h3_block(M,'ms_ext_primary')
res['H4_ext'],_,_,_=H4mod.evaluate(M,'ms_ext_primary',nperm=200)
Fall=H4mod.features(M); tr=M[M.ms_ext_primary.notna()&(M.train==1)]; ytr=tr.ms_ext_primary.astype(int).values
mods={nm:H4mod.fit_pred(Fall.loc[tr.index],ytr,Fall.loc[tr.index[:2]],cols)[1:] for nm,cols in [('full',H4mod.FULL),('B0',H4mod.B0),('B1',H4mod.B1)]}
cols={'full':H4mod.FULL,'B0':H4mod.B0,'B1':H4mod.B1}
def score(nm,idx):
    lr,(med,mu,sd)=mods[nm]; X=Fall.loc[idx,cols[nm]].fillna(med); return lr.predict_proba((X-mu)/sd)[:,1]
for layer in ['gs','hb']:
    lab=f'ms_{layer}_primary'
    res[f'H3_{layer}'],_=H3mod.h3_block(M,lab)
    te=M[M[lab].notna()&(M.train==0)]; y=te[lab].astype(int).values
    pf,p0,p1=score('full',te.index),score('B0',te.index),score('B1',te.index)
    res[f'H4_{layer}']=dict(n_test=int(len(te)),pos=int(y.sum()),vsB0=H4mod.delong(y,pf,p0),vsB1=H4mod.delong(y,pf,p1))
res['A1_train_n']=int(len(tr)); res['A1_coef_full']=dict(zip(H4mod.FULL,map(float,mods['full'][0].coef_[0])))
dump(res,'ph_ms.json'); print(json.dumps({k:res[k] for k in res if k.startswith('n_')}))
