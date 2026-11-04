"""Shared helpers for C10+ analyses (definitions per frozen pre-registration v2.0, SHA256 126f515f...)."""
import numpy as np, pandas as pd, statsmodels.api as sm, json, warnings
from scipy.stats import norm, fisher_exact
warnings.filterwarnings('ignore')
W="./work/"
R="./results/"
SEED=20261103
HLA=('6',28510120,33480577)
def load():
    # NOTE: OT eaf is missing for UKB-PPP; MAF = gnomAD NFE AF of the lead variant (documented deviation)
    M=pd.read_parquet(W+"master_cs.parquet")
    lp=M.log10p.astype(float); lp=np.where(lp<0,-lp,lp)   # OT log10p sign convention guard
    M['mlog10p']=lp; M['mlp']=np.minimum(lp,300)/10
    M['absb']=M.beta.abs(); M['maf']=np.minimum(M.lead_nfe_af,1-M.lead_nfe_af); M['logn95']=np.log(np.maximum(M.n95.fillna(1),1))
    M['dist_tss']=(M.pos_b38-M.tss).abs(); M['logdist']=np.log10(1+M.dist_tss.fillna(1e6))
    M['PAV']=M.pav.fillna(False).astype(int); M['PAV_canon']=M.pav_canon.fillna(False).astype(int)
    M['PAV_pp02']=(M.pav_pp.fillna(0)>=0.2).astype(int); M['PAV_lead']=M.pav_lead.fillna(False).astype(int)
    M['v3']=M.v3.fillna(False).astype(bool); M['v4']=M.v4.fillna(False).astype(bool)
    M['chrom']=M.chrom.astype(str)
    M['hla']=((M.chrom==HLA[0])&(M.pos_b38>=HLA[1])&(M.pos_b38<=HLA[2])).astype(int)
    M['hotspot']=(M.n_trans_genes_shared>=5).astype(int)
    M['train']=M.chrom.apply(lambda c: (c=='X') or (c.isdigit() and int(c)%2==1)).astype(int)
    return M
def derive(M,thr=0.8,var=''):
    """var: '' primary prior; '_p12_5e6'; '_p12_1e6'; '_inter'"""
    M=M.copy()
    g=lambda c: M[f'{c}_h4{var}']
    M['A_gtex']=(g('gtex')>=thr).astype(int); M['A_nongtex']=(g('nongtex')>=thr).astype(int)
    M['A_eqtl']=((M.A_gtex+M.A_nongtex)>0).astype(int); M['A_iso']=(g('iso')>=thr).astype(int)
    M['W_gtex']=(g('wrong')>=thr).astype(int); M['rep_sun']=(g('sun')>=thr).astype(int)
    M['concordant']=((M.v3&(M.rep_sun==1))|(M.v4&(M.rep_decode==1))).astype(int)
    M['conc_cat']=np.where(M.concordant==1,'concordant',np.where(M.v3|M.v4,'discordant','na'))
    sd_sun=np.sign(M.sun_b_ukb)*np.sign(M.sun_b_sun); sd_dec=np.sign(M.dec_b_ukb)*np.sign(M.dec_b_dec)
    opp=((M.v3&(M.rep_sun==1)&(sd_sun<0))|(M.v4&(M.rep_decode==1)&(sd_dec<0)))
    samedir_ok=((M.v3&(M.rep_sun==1)&(sd_sun>0))|(M.v4&(M.rep_decode==1)&(sd_dec>0)))
    M['concordant_dir']=(samedir_ok|((M.concordant==1)&~opp&~samedir_ok)).astype(int)   # opposite-direction pairs not counted as concordant
    return M
def glm(df,y,xs,cluster='ga',terms=None):
    d=df.dropna(subset=[y,*xs]).copy()
    X=sm.add_constant(d[xs].astype(float),has_constant='add'); yy=d[y].astype(float)
    fam=sm.families.Binomial()
    m0=sm.GLM(yy,X,family=fam).fit()
    try: mc=sm.GLM(yy,X,family=fam).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(d[cluster])[0]})
    except Exception: mc=m0
    out={}
    for t in (terms or xs[:1]):
        b=mc.params[t]; se=mc.bse[t]; se0=m0.bse[t]
        out[t]=dict(OR=float(np.exp(b)),lo=float(np.exp(b-1.96*se)),hi=float(np.exp(b+1.96*se)),P=float(mc.pvalues[t]),P_model=float(m0.pvalues[t]),
                    lo_model=float(np.exp(b-1.96*se0)),hi_model=float(np.exp(b+1.96*se0)),n=int(len(d)),n_pos=int(yy.sum()),n_exposed=int(d[t].sum()) if set(d[t].unique())<= {0,1} else None,n_genes=int(d[cluster].nunique()))
    return out,mc
def diff_test(mc,a,b):
    dl=mc.params[a]-mc.params[b]; V=mc.cov_params(); sed=np.sqrt(V.loc[a,a]+V.loc[b,b]-2*V.loc[a,b])
    return dict(diff_logOR=float(dl),se=float(sed),P_one_sided=float(norm.sf(dl/sed)),OR_a=float(np.exp(mc.params[a])),OR_b=float(np.exp(mc.params[b])))
def fisher(df,g,y):
    t=pd.crosstab(df[g],df[y]).reindex(index=[1,0],columns=[1,0],fill_value=0)
    OR,p=fisher_exact(t.values)
    return dict(table=t.values.tolist(),OR=float(OR),P=float(p),rate_g1=float(t.loc[1,1]/max(t.loc[1].sum(),1)),rate_g0=float(t.loc[0,1]/max(t.loc[0].sum(),1)),n_g1=int(t.loc[1].sum()),n_g0=int(t.loc[0].sum()))
def pav_match(df,y,strat,covs=('mlp','absb','maf','logn95'),k=3,caliper=0.5,pav='PAV',seed=SEED):
    """1:k nearest-neighbour Mahalanobis matching without replacement within strata; returns matched df with set ids."""
    rng=np.random.default_rng(seed); d=df.dropna(subset=list(covs)+[y]).copy()
    Z=(d[list(covs)]-d[list(covs)].mean())/d[list(covs)].std(); VI=np.linalg.pinv(np.cov(Z.values.T))
    d[[c+'_z' for c in covs]]=Z.values
    sets=[]; sid=0
    for s,ds in d.groupby(strat):
        tr=ds[ds[pav]==1]; co=ds[ds[pav]==0].copy()
        if len(tr)==0 or len(co)==0: continue
        Zc=co[[c+'_z' for c in covs]].values; avail=np.ones(len(co),bool)
        for i in rng.permutation(len(tr)):
            zt=tr[[c+'_z' for c in covs]].values[i]
            dd=Zc-zt; dist=np.sqrt(np.einsum('ij,jk,ik->i',dd,VI,dd)); dist[~avail]=np.inf
            idx=np.argsort(dist)[:k]; idx=idx[dist[idx]<=caliper]
            if len(idx)==0: continue
            avail[idx]=False
            sets.append(tr.iloc[[i]].assign(set_id=sid)); sets.append(co.iloc[idx].assign(set_id=sid)); sid+=1
    if not sets: return None
    Mt=pd.concat(sets)
    return Mt
def smd(df,g,covs=('mlp','absb','maf','logn95')):
    out={}
    for c in covs:
        a=df[df[g]==1][c]; b=df[df[g]==0][c]; out[c]=float((a.mean()-b.mean())/np.sqrt((a.var()+b.var())/2))
    return out
def fisher_z_diff(r1,n1,r2,n2):
    z=(np.arctanh(r1)-np.arctanh(r2))/np.sqrt(1/(n1-3)+1/(n2-3)); return float(z), float(norm.sf(-z))  # one-sided P that r1<r2
def holm(pdict):
    items=sorted(pdict.items(),key=lambda x:x[1]); m=len(items); out={}; run=0
    for i,(k,p) in enumerate(items):
        adj=min(1,(m-i)*p); run=max(run,adj); out[k]=run
    return out
def dump(obj,name):
    json.dump(obj,open(R+name,'w'),indent=1,default=lambda o: o.item() if hasattr(o,'item') else str(o))
