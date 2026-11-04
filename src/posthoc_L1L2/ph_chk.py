import sys,os,numpy as np,pandas as pd
PH='./posthoc_L1L2'; sys.path.insert(0,PH); os.environ['PH_TAG']='L2'
from c10common import *
import h4 as H4mod
from sklearn.metrics import roc_auc_score
A4='./a4/'
M=derive(load()).merge(pd.read_parquet(A4+'a4_cs_fields.parquet')[['a','ms_gs_primary','ms_hb_primary']],on='a',how='left').reset_index(drop=True)
F=H4mod.features(M)
for lab in ['ms_gs_primary','ms_hb_primary']:
    T=M[M[lab].notna()]; y=T[lab].astype(int).values
    p,lr,_=H4mod.fit_pred(F.loc[T.index],y,F.loc[T.index],H4mod.B0)
    print(lab,'in-sample B0 AUC',round(roc_auc_score(y,p),3),'coef',lr.coef_[0].round(3),pd.crosstab(T.conc_cat,y).to_dict())
    # per-chromosome concordance association heterogeneity
    r=T.assign(y=y).groupby(['chrom','conc_cat']).y.mean().unstack(); print(r.round(2).head(25).to_string())
