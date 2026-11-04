"""Addendum A1: Seer full summary statistics -> MS support labels (ms_ext_*) and H3i-ext / H3ii-ext / H4-ext.
Pre-registration addendum: C10plus_prereg_addendum_A1_seer_full_frozen.md (SHA256 6d54bd54...)."""
import sys, glob, os, json, numpy as np, pandas as pd, duckdb
sys.path.insert(0,'.')
from c10common import *
from scipy.stats import norm
import h3 as H3mod, h4 as H4mod
SE="./seer_ext/"
B="../data_v1/db/parquet/"
mode=sys.argv[1] if len(sys.argv)>1 else 'full'
ST=pd.read_csv(SE+'seer_study_windows.csv',dtype={'chrom':str})
L=pd.read_parquet(W+'ukb_cis_locus_grch37.parquet'); L=L[L.is95&L.pos37.notna()].copy().rename(columns={'beta':'b_ukb_locus'}); L['pos37']=L.pos37.astype(int); L['chrom']=L.chrom.astype(str)
L['a1']=np.minimum(L.ref,L.alt); L['a2']=np.maximum(L.ref,L.alt)
files=sorted(glob.glob(SE+'all/*.cis.tsv.gz'))
def readf(f):
    d=pd.read_csv(f,sep='\t',dtype={'chromosome':str},usecols=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','n'])
    d['gcst']=os.path.basename(f).split('.')[0]; return d
if mode=='offset':
    # offset determination using positions/alleles only (no association values inspected)
    res={}
    for off in (0,1):
        hit=tot=0
        for f in files[:40]:
            d=pd.read_csv(f,sep='\t',dtype={'chromosome':str},usecols=['chromosome','base_pair_location','effect_allele','other_allele'])
            d['pos']=d.base_pair_location+off; d['a1']=np.minimum(d.effect_allele,d.other_allele); d['a2']=np.maximum(d.effect_allele,d.other_allele)
            g=ST.set_index('gcst').loc[os.path.basename(f).split('.')[0],'gene_id']
            Lg=L[L.ga==g]
            m=Lg.merge(d,left_on=['chrom','pos37','a1','a2'],right_on=['chromosome','pos','a1','a2'])
            hit+=m.variant_id.nunique(); tot+=Lg.variant_id.nunique()
        res[off]=dict(matched=int(hit),ukb_cs_variants=int(tot),rate=hit/max(tot,1))
    print(json.dumps(res)); json.dump(res,open(SE+'offset_check.json','w')); sys.exit()
OFF=json.load(open(SE+'offset_choice.json'))['offset']
con=duckdb.connect()
cs=con.sql(f"select studyLocusId a, nSamples, sampleSize from '{B}qtl_cs.parquet' where projectId='UKB_PPP_EUR' and not isTransQtl").df()
M=derive(load()); M=M.merge(cs,on='a',how='left'); M['N_ukb']=M.nSamples.fillna(M.sampleSize)
M['z_ukb']=norm.isf(np.minimum(10.0**(-np.minimum(M.mlog10p,300)),1.0)/2)
pairs=[]
st_gene=ST.set_index('gcst').gene_id.to_dict()
for f in files:
    d=readf(f)
    if len(d)==0: continue
    g=st_gene.get(d.gcst.iloc[0])
    Lg=L[L.ga==g]
    if len(Lg)==0: continue
    d['pos']=d.base_pair_location+OFF; d['a1']=np.minimum(d.effect_allele,d.other_allele); d['a2']=np.maximum(d.effect_allele,d.other_allele)
    m=Lg.merge(d,left_on=['chrom','pos37','a1','a2'],right_on=['chromosome','pos','a1','a2'])
    if len(m)==0: continue
    maxn=d.n.max()
    m['maf_s']=np.minimum(m.effect_allele_frequency,1-m.effect_allele_frequency)
    pal=m.apply(lambda r: {r.ref,r.alt} in ({'A','T'},{'C','G'}),axis=1)
    m=m[(m.maf_s>=0.01)&(m.n>=0.8*maxn)&~(pal&(m.effect_allele_frequency.between(0.4,0.6)))].copy()
    if len(m)==0: continue
    m['sgn']=np.where(m.effect_allele==m.alt,1,np.where(m.effect_allele==m.ref,-1,0)); m=m[m.sgn!=0]
    m['b_s']=m.sgn*m.beta
    pairs.append(m.assign(gene=g))
P=pd.concat(pairs)
P=P.merge(M[['a','lead_variant_id','beta','z_ukb','N_ukb']].rename(columns={'beta':'lead_beta'}),on='a')
P['b_ukb']=np.where(P.variant_id==P.lead_variant_id,P.lead_beta,P.b_ukb_locus)
# index variant per (cs, study)
P=P.sort_values(['a','gcst','pp'],ascending=[True,True,False])
n_eval=P.groupby(['a','gcst']).size().rename('n_eval')
I=P.drop_duplicates(['a','gcst']).merge(n_eval,left_on=['a','gcst'],right_index=True)
K=len(I); alpha=0.05/K; zc=norm.isf(alpha/2); zc05=norm.isf(0.025)
I['z_exp']=I.z_ukb*np.sqrt(I.n/I.N_ukb)
I['power']=norm.cdf(I.z_exp-zc)+norm.cdf(-I.z_exp-zc); I['power05']=norm.cdf(I.z_exp-zc05)+norm.cdf(-I.z_exp-zc05)
same=np.sign(I.b_s)==np.sign(I.b_ukb)
I['sup']=(I.p_value<alpha)&same; I['opp_sig']=(I.p_value<alpha)&~same
I['lab']=np.where(I.sup,1.0,np.where(I.power>=0.99,0.0,np.nan))
I['sup05']=(I.p_value<0.05)&same; I['lab05']=np.where(I.sup05,1.0,np.where(I.power05>=0.99,0.0,np.nan))
# cs-wide label
P2=P.merge(n_eval,left_on=['a','gcst'],right_index=True); P2['sw']=(P2.p_value<alpha/P2.n_eval)&(np.sign(P2.b_s)==np.sign(P2.b_ukb))
cw=P2.groupby(['a','gcst']).sw.any().rename('sup_cw'); I=I.merge(cw,left_on=['a','gcst'],right_index=True)
I['lab_cw']=np.where(I.sup_cw,1.0,np.where(I.power>=0.99,0.0,np.nan))
def comb(s):
    s=s.dropna(); 
    if (s==1).any(): return 1.0
    if (s==0).any(): return 0.0
    return np.nan
agg=I.groupby('a').agg(ms_ext_primary=('lab',comb),ms_ext_nominal=('lab05',comb),ms_ext_cswide=('lab_cw',comb),seer_ext_opposite_sig=('opp_sig','any'),seer_ext_studies=('gcst',lambda s:';'.join(sorted(s))))
best=I.sort_values(['a','power'],ascending=[True,False]).drop_duplicates('a').set_index('a')[['variant_id','b_s','standard_error','p_value','power','z_exp']]
best.columns=['seer_ext_index_variant','seer_ext_beta_aligned','seer_ext_se','seer_ext_p','seer_ext_power','seer_ext_zexp']
E=agg.join(best).reset_index()
E.to_parquet(SE+'ms_ext_labels.parquet',index=False); I.to_parquet(SE+'ms_ext_pairs.parquet',index=False)
info=dict(offset=OFF,n_files=len(files),K_pairs=int(K),alpha=float(alpha),n_cs_assessed=int(len(E)),labels=E.ms_ext_primary.value_counts(dropna=False).to_dict(),
          nominal=E.ms_ext_nominal.value_counts(dropna=False).to_dict(),cswide=E.ms_ext_cswide.value_counts(dropna=False).to_dict(),opp_sig=int(E.seer_ext_opposite_sig.sum()))
# ---------------- endpoints ----------------
M=M.merge(E,on='a',how='left')
res={'info':info}
res['H3_ext'],T=H3mod.h3_block(M,'ms_ext_primary')
res['H3_ext_nominal'],_=H3mod.h3_block(M,'ms_ext_nominal',allow_stop=False)
res['H3_ext_cswide'],_=H3mod.h3_block(M,'ms_ext_cswide',allow_stop=False)
new=M[M.ms_sheet.isna()]
res['H3_ext_new'],_=H3mod.h3_block(new,'ms_ext_primary',allow_stop=False)
lp=M.log10p.abs(); z=norm.isf(np.minimum(10**(-lp.clip(upper=300))/2,0.5)); z=np.where(lp>=300,40,z)
ratio=M.beta.abs()/(z/np.sqrt(2*M.maf*(1-M.maf)*34000)); M['flag_beta_inconsistent']=(ratio>5).astype(int)
res['H3_ext_excl_beta_flag'],_=H3mod.h3_block(M[M.flag_beta_inconsistent==0],'ms_ext_primary',allow_stop=False)
res['H4_ext'],Tm,Fm,pm=H4mod.evaluate(M,'ms_ext_primary')
# H4-ext-new: same model (trained on all ext training rows) evaluated on ext-new held-out rows
Fall=H4mod.features(M); tr=M[M.ms_ext_primary.notna()&(M.train==1)]; te=M[M.ms_ext_primary.notna()&(M.train==0)&M.ms_sheet.isna()]
pr=[H4mod.fit_pred(Fall.loc[tr.index],tr.ms_ext_primary.astype(int).values,Fall.loc[te.index],cols)[0] for cols in (H4mod.FULL,H4mod.B0)]
res['H4_ext_new']=dict(n_test=int(len(te)),pos_test=int(te.ms_ext_primary.sum()),**H4mod.delong(te.ms_ext_primary.astype(int).values,pr[0],pr[1]))
# agreement with MSPA primary
both=M[M.ms_primary.notna()&M.ms_ext_primary.notna()]
if len(both):
    from sklearn.metrics import cohen_kappa_score
    res['agreement_mspa_vs_ext']=dict(n=int(len(both)),agree=float((both.ms_primary==both.ms_ext_primary).mean()),kappa=float(cohen_kappa_score(both.ms_primary,both.ms_ext_primary)),table=pd.crosstab(both.ms_primary,both.ms_ext_primary).values.tolist())
M['tier']=H4mod.tiers(M); d=M[M.ms_ext_primary.notna()&(M.tier>0)]
res['tiers_ext']=dict(rates=d.groupby('tier').ms_ext_primary.agg(['sum','size','mean']).reset_index().to_dict('records'),trend=H4mod.cochran_armitage(d.tier,d.ms_ext_primary))
res['holm_addendum']=holm({'H3i-ext':res['H3_ext']['H3i_logit']['A_eqtl']['P'],'H3ii-ext':res['H3_ext']['H3ii_CNP_vs_CA']['P'],'H4-ext':res['H4_ext']['delong_full_vs_B0']['P_one_sided']})
T.to_csv(SE+'h3_ext_table.csv',index=False)
pd.DataFrame({'y':Tm[Tm.train==0].ms_ext_primary.values,'p_full':pm['full'],'p_B0':pm['B0'],'p_B1':pm['B1']}).to_csv(SE+'h4_ext_heldout_preds.csv',index=False)
dump(res,'h3h4_ext.json'); print(json.dumps(info,default=str))
