"""H4: robustness grade trained on odd chromosomes (+X), evaluated on held-out even chromosomes."""
import sys; sys.path.insert(0,'.')
from c10common import *
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2
def delong(y,p1,p2):
    """Fast DeLong (Sun & Xu 2014) for two correlated AUCs; returns auc1, auc2, z, one-sided P(auc1>auc2)."""
    y=np.asarray(y); pos=np.where(y==1)[0]; neg=np.where(y==0)[0]; m,n=len(pos),len(neg)
    def comps(p):
        V10=np.array([((p[i]>p[neg]).mean()+0.5*(p[i]==p[neg]).mean()) for i in pos])
        V01=np.array([((p[pos]>p[j]).mean()+0.5*(p[pos]==p[j]).mean()) for j in neg])
        return V10.mean(),V10,V01
    a1,v10a,v01a=comps(np.asarray(p1)); a2,v10b,v01b=comps(np.asarray(p2))
    S10=np.cov(np.vstack([v10a,v10b])); S01=np.cov(np.vstack([v01a,v01b])); S=S10/m+S01/n
    var=S[0,0]+S[1,1]-2*S[0,1]; z=(a1-a2)/np.sqrt(var) if var>0 else np.nan
    return dict(auc1=float(a1),auc2=float(a2),dauc=float(a1-a2),z=float(z),P_one_sided=float(norm.sf(z)),se_auc1=float(np.sqrt(S[0,0])),se_auc2=float(np.sqrt(S[1,1])))
def features(M):
    F=pd.DataFrame(index=M.index)
    F['A_eqtl']=M.A_eqtl; F['log1p_nanchor']=np.log1p(M.n_anchor_datasets); F['A_iso']=M.A_iso; F['PAV']=M.PAV; F['pav_pp']=M.pav_pp.fillna(0)
    F['PAVxNoAnchor']=M.PAV*(1-M.A_eqtl)
    F['conc_concordant']=(M.conc_cat=='concordant').astype(int); F['conc_discordant']=(M.conc_cat=='discordant').astype(int)
    F['mlp']=M.mlp; F['absb']=M.absb; F['maf']=M.maf; F['logn95']=M.logn95; F['purityMinR2']=M.purityMinR2; F['logdist']=M.logdist
    return F
FULL=['A_eqtl','log1p_nanchor','A_iso','PAV','pav_pp','PAVxNoAnchor','conc_concordant','conc_discordant','mlp','absb','maf','logn95','purityMinR2','logdist']
B0=['conc_concordant','conc_discordant']; B1=B0+['mlp','absb','maf']
def fit_pred(Ftr,ytr,Fte,cols):
    med=Ftr[cols].median(); Xtr=Ftr[cols].fillna(med); Xte=Fte[cols].fillna(med)
    mu=Xtr.mean(); sd=Xtr.std().replace(0,1)
    lr=LogisticRegression(C=1.0,penalty='l2',max_iter=5000).fit((Xtr-mu)/sd,ytr)
    return lr.predict_proba((Xte-mu)/sd)[:,1], lr, (med,mu,sd)
def calib(y,p):
    p=np.clip(p,1e-6,1-1e-6); lp=np.log(p/(1-p))
    m=sm.GLM(y,sm.add_constant(lp),family=sm.families.Binomial()).fit()
    mi=sm.GLM(y,np.ones((len(y),1)),family=sm.families.Binomial(),offset=lp).fit()
    return dict(slope=float(m.params[1]),intercept_citl=float(mi.params[0]),brier=float(brier_score_loss(y,p)))
def evaluate(M,target,nperm=1000):
    T=M[M[target].notna()].copy(); F=features(T); y=T[target].astype(int).values
    tr=T.train.values==1; te=~tr
    o=dict(target=target,n_train=int(tr.sum()),n_test=int(te.sum()),pos_train=int(y[tr].sum()),pos_test=int(y[te].sum()))
    if te.sum()<40 or min(y[te].sum(),(1-y[te]).sum())<10: o['stop_rule_triggered']=True
    else: o['stop_rule_triggered']=False
    preds={}
    for nm,cols in [('full',FULL),('B0',B0),('B1',B1)]:
        preds[nm],lr,_=fit_pred(F[tr],y[tr],F[te],cols)
        o[f'auc_{nm}']=float(roc_auc_score(y[te],preds[nm])); o[f'calib_{nm}']=calib(y[te],preds[nm])
        if nm=='full': o['coef_full_std']=dict(zip(cols,map(float,lr.coef_[0])))
    o['delong_full_vs_B0']=delong(y[te],preds['full'],preds['B0']); o['delong_full_vs_B1']=delong(y[te],preds['full'],preds['B1'])
    rng=np.random.default_rng(SEED); null=[]
    for i in range(nperm):
        yp=rng.permutation(y[tr]); pp,_,_=fit_pred(F[tr],yp,F[te],FULL); null.append(roc_auc_score(y[te],pp))
    null=np.array(null); o['perm_null_auc_mean']=float(null.mean()); o['perm_null_auc_q95']=float(np.quantile(null,0.95)); o['perm_P']=float((1+(null>=o['auc_full']).sum())/(1+nperm))
    return o,T,F,preds
def tiers(M):
    t=np.select([(M.A_eqtl==1)&(M.PAV==0),((M.A_eqtl==1)&(M.PAV==1))|((M.A_eqtl==0)&(M.A_iso==1)&(M.PAV==0)),(M.A_eqtl==0)&(M.A_iso==0)&(M.PAV==0),(M.A_eqtl==0)&(M.PAV==1)],[1,2,3,4],0)
    return t
def cochran_armitage(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(y); pbar=y.mean()
    T=np.sum(x*(y-pbar)); V=pbar*(1-pbar)*np.sum((x-x.mean())**2)
    z=T/np.sqrt(V); return dict(z=float(z),P_two_sided=float(2*norm.sf(abs(z))))
if __name__=="__main__":
    M=derive(load()); M['tier']=tiers(M); res={}
    res['H4_ms'],Tm,Fm,pm=evaluate(M,'ms_primary')
    res['H4_xa'],Tx,Fx,px=evaluate(M,'xa_rep',nperm=200)
    # transfer: MS-trained models scored on held-out xa_rep
    Fall=features(M); Tms=M[M.ms_primary.notna()&(M.train==1)]; yms=Tms.ms_primary.astype(int).values
    te=M[(M.train==0)&M.xa_rep.notna()]
    tr_ms=[fit_pred(Fall.loc[Tms.index],yms,Fall.loc[te.index],cols)[0] for cols in (FULL,B0)]
    res['transfer_ms_to_xa']=delong(te.xa_rep.astype(int).values,tr_ms[0],tr_ms[1])
    # also MS-trained full model on held-out MS secondary labels
    for lab in ['ms_mspa','ms_seer']:
        tt=M[(M.train==0)&M[lab].notna()]
        pr=[fit_pred(Fall.loc[Tms.index],yms,Fall.loc[tt.index],cols)[0] for cols in (FULL,B0)]
        res[f'transfer_ms_to_{lab}']=delong(tt[lab].astype(int).values,pr[0],pr[1])
    # rule-based tiers: all chromosomes
    tr={}
    for lab in ['ms_primary','ms_mspa','ms_seer','xa_rep']:
        d=M[M[lab].notna()&(M.tier>0)]
        tr[lab]=dict(rates=d.groupby('tier')[lab].agg(['sum','size','mean']).reset_index().to_dict('records'),trend=cochran_armitage(d.tier,d[lab]))
    res['tiers_vs_truth']=tr; res['tier_counts']=M.tier.value_counts().sort_index().to_dict()
    # final model refit on all chromosomes (MS target) -> confidence for every cs
    Tall=M[M.ms_primary.notna()]; p_all,lr,_=fit_pred(Fall.loc[Tall.index],Tall.ms_primary.astype(int).values,Fall,FULL)
    res['final_model_coef_std']=dict(zip(FULL,map(float,lr.coef_[0])))
    Tx_all=M[M.xa_rep.notna()]; p_xa,lr2,_=fit_pred(Fall.loc[Tx_all.index],Tx_all.xa_rep.astype(int).values,Fall,FULL)
    res['final_xa_model_coef_std']=dict(zip(FULL,map(float,lr2.coef_[0])))
    out=M[['a']].copy(); out['tier']=M.tier; out['p_ms_support_model']=p_all; out['p_xa_rep_model']=p_xa
    out.to_parquet(R+'h4_scores.parquet',index=False)
    pd.DataFrame({'y':Tm[Tm.train==0].ms_primary.values,'p_full':pm['full'],'p_B0':pm['B0'],'p_B1':pm['B1']}).to_csv(R+'h4_heldout_ms_preds.csv',index=False)
    pd.DataFrame({'y':Tx[Tx.train==0].xa_rep.values,'p_full':px['full'],'p_B0':px['B0'],'p_B1':px['B1']}).to_csv(R+'h4_heldout_xa_preds.csv',index=False)
    dump(res,'h4.json')
    print(json.dumps({k:res['H4_ms'][k] for k in ['n_train','n_test','pos_test','auc_full','auc_B0','auc_B1','delong_full_vs_B0','perm_P','stop_rule_triggered']},indent=1))
    print(json.dumps({k:res['H4_xa'][k] for k in ['n_train','n_test','auc_full','auc_B0','auc_B1','delong_full_vs_B0']},indent=1)); print(json.dumps(res['tiers_vs_truth'],indent=0,default=str)[:3000])
