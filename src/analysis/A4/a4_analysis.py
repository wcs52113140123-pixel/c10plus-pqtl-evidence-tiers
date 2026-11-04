"""Addendum A4: independent MS truth layers (GS primary, HOLBAEK secondary). Prereg A4 v1.0 SHA256 cc2da044...fd7.
Label rule identical to A1 (ms_ext.py); frozen A1 model for H4-A4."""
import sys, glob, os, json, numpy as np, pandas as pd, duckdb
sys.path.insert(0,'.')
from c10common import *
from scipy.stats import norm, fisher_exact
from sklearn.metrics import cohen_kappa_score
import h3 as H3mod, h4 as H4mod
A4='./a4/'; SE='./seer_ext/'
B="../data_v1/db/parquet/"
mode=sys.argv[1] if len(sys.argv)>1 else 'full'
L0=pd.read_parquet(W+'ukb_cis_locus_grch37.parquet'); L0=L0[L0.is95].copy().rename(columns={'beta':'b_ukb_locus'}); L0['chrom']=L0.chrom.astype(str)
L0['a1']=np.minimum(L0.ref,L0.alt); L0['a2']=np.maximum(L0.ref,L0.alt)
GSmap=pd.read_csv(A4+'gs_entry_gene_windows.tsv',sep='\t',dtype={'chrom':str})
HBmap=pd.read_csv(A4+'hb_jobs.tsv',sep='\t',header=None,names=['gcst','label','gene_id','chrom','w0','w1','url'],dtype={'chrom':str})
def read_layer(layer,f):
    if layer=='gs':
        d=pd.read_csv(f,sep='\t',dtype={'CHR':str})
        d=d.rename(columns={'CHR':'chromosome','POS':'base_pair_location','A1':'effect_allele','A2':'other_allele','AF1':'effect_allele_frequency','BETA':'beta','SE':'standard_error','P':'p_value','N':'n'})
        d['study']=os.path.basename(f).replace('.cis.tsv.gz','')+'.fastGWA.gz'
    else:
        d=pd.read_csv(f,sep='\t',dtype={'chromosome':str})
        d.columns=[c.strip() for c in d.columns]
        d['study']=os.path.basename(f).split('.')[0]
        d=d.drop(columns=[c_ for c_ in ['variant_id','rsid','ref_allele','hm_coordinate_conversion','hm_code','trait(UniprotID_GeneName)'] if c_ in d.columns])
        if 'n' not in d.columns or d.n.isna().all(): d['n']=1909
        d['n']=pd.to_numeric(d.n,errors='coerce').fillna(1909)
    for c_ in ['beta','standard_error','p_value','effect_allele_frequency','base_pair_location']: d[c_]=pd.to_numeric(d[c_],errors='coerce')
    d=d.dropna(subset=['base_pair_location','beta','p_value'])
    d['effect_allele']=d.effect_allele.astype(str).str.upper(); d['other_allele']=d.other_allele.astype(str).str.upper()
    d['a1']=np.minimum(d.effect_allele,d.other_allele); d['a2']=np.maximum(d.effect_allele,d.other_allele)
    return d
def study_genes(layer,multi=False):
    if layer=='gs':
        g=GSmap if multi else GSmap[GSmap.n_genes==1]
        return g.groupby('name').gene_id.apply(list).to_dict()
    return {r.gcst:[r.gene_id] for r in HBmap.itertuples()}
FILES={'gs':sorted(glob.glob(A4+'gs_cis/*.cis.tsv.gz')),'hb':sorted(glob.glob(A4+'hb_cis/*.cis.tsv.gz'))}
POS={'gs':'pos37','hb':'pos_b38'}
def offset_check(layer):
    sg=study_genes(layer,True); res={}
    for off in (0,1):
        hit=tot=0
        for f in FILES[layer][:40]:
            d=read_layer(layer,f)
            if len(d)==0: continue
            d=d[['chromosome','base_pair_location','a1','a2','study']]
            d['pos']=d.base_pair_location+off
            for g in sg.get(d.study.iloc[0],[]):
                Lg=L0[L0.ga==g].dropna(subset=[POS[layer]]).copy(); Lg['pos']=Lg[POS[layer]].astype(int)
                m=Lg.merge(d,left_on=['chrom','pos','a1','a2'],right_on=['chromosome','pos','a1','a2'])
                hit+=m.variant_id.nunique(); tot+=Lg.variant_id.nunique()
        res[off]=dict(matched=int(hit),ukb_cs_variants=int(tot),rate=hit/max(tot,1))
    return res
if mode=='offset':
    out={l:offset_check(l) for l in ['gs','hb']}; json.dump(out,open(A4+'offset_check.json','w'),indent=1); print(json.dumps(out)); sys.exit()
OFFC=json.load(open(A4+'offset_check.json')); OFF={l:max((0,1),key=lambda o:OFFC[l][str(o)]['rate']) for l in OFFC}
con=duckdb.connect()
cs=con.sql(f"select studyLocusId a, nSamples, sampleSize from '{B}qtl_cs.parquet' where projectId='UKB_PPP_EUR' and not isTransQtl").df()
M=derive(load()); M=M.merge(cs,on='a',how='left'); M['N_ukb']=M.nSamples.fillna(M.sampleSize)
M['z_ukb']=norm.isf(np.minimum(10.0**(-np.minimum(M.mlog10p,300)),1.0)/2)
Ext=pd.read_parquet(SE+'ms_ext_labels.parquet'); M=M.merge(Ext[['a','ms_ext_primary','seer_ext_studies']],on='a',how='left')
def labels(layer,multi=False):
    sg=study_genes(layer,multi); pairs=[]
    for f in FILES[layer]:
        d=read_layer(layer,f)
        if len(d)==0: continue
        st=d.study.iloc[0]
        for g in sg.get(st,[]):
            Lg=L0[L0.ga==g].dropna(subset=[POS[layer]]).copy()
            if len(Lg)==0: continue
            Lg['pos']=Lg[POS[layer]].astype(int); dd=d.copy(); dd['pos']=dd.base_pair_location.astype(int)+OFF[layer]
            m=Lg.merge(dd,left_on=['chrom','pos','a1','a2'],right_on=['chromosome','pos','a1','a2'])
            if len(m)==0: continue
            maxn=d.n.max(); m['maf_s']=np.minimum(m.effect_allele_frequency,1-m.effect_allele_frequency)
            pal=m.apply(lambda r: {r.ref,r.alt} in ({'A','T'},{'C','G'}),axis=1)
            m=m[(m.maf_s>=0.01)&(m.n>=0.8*maxn)&~(pal&(m.effect_allele_frequency.between(0.4,0.6)))].copy()
            if len(m)==0: continue
            m['sgn']=np.where(m.effect_allele==m.alt,1,np.where(m.effect_allele==m.ref,-1,0)); m=m[m.sgn!=0]
            m['b_s']=m.sgn*m.beta; pairs.append(m.assign(gene=g,study=st))
    P=pd.concat(pairs)
    P=P.merge(M[['a','lead_variant_id','beta','z_ukb','N_ukb']].rename(columns={'beta':'lead_beta'}),on='a')
    P['b_ukb']=np.where(P.variant_id==P.lead_variant_id,P.lead_beta,P.b_ukb_locus)
    P=P.sort_values(['a','study','pp'],ascending=[True,True,False])
    n_eval=P.groupby(['a','study']).size().rename('n_eval')
    I=P.drop_duplicates(['a','study']).merge(n_eval,left_on=['a','study'],right_index=True)
    K=len(I); alpha=0.05/K; zc=norm.isf(alpha/2); zc05=norm.isf(0.025)
    I['z_exp']=I.z_ukb*np.sqrt(I.n/I.N_ukb)
    I['power']=norm.cdf(I.z_exp-zc)+norm.cdf(-I.z_exp-zc); I['power05']=norm.cdf(I.z_exp-zc05)+norm.cdf(-I.z_exp-zc05)
    same=np.sign(I.b_s)==np.sign(I.b_ukb)
    I['sup']=(I.p_value<alpha)&same; I['opp_sig']=(I.p_value<alpha)&~same
    I['lab']=np.where(I.sup,1.0,np.where(I.power>=0.99,0.0,np.nan))
    I['sup05']=(I.p_value<0.05)&same; I['lab05']=np.where(I.sup05,1.0,np.where(I.power05>=0.99,0.0,np.nan))
    P2=P.merge(n_eval,left_on=['a','study'],right_index=True); P2['sw']=(P2.p_value<alpha/P2.n_eval)&(np.sign(P2.b_s)==np.sign(P2.b_ukb))
    cw=P2.groupby(['a','study']).sw.any().rename('sup_cw'); I=I.merge(cw,left_on=['a','study'],right_index=True)
    I['lab_cw']=np.where(I.sup_cw,1.0,np.where(I.power>=0.99,0.0,np.nan))
    def comb(s):
        s=s.dropna()
        if (s==1).any(): return 1.0
        if (s==0).any(): return 0.0
        return np.nan
    p=f'{layer}_'
    agg=I.groupby('a').agg(**{f'ms_{layer}_primary':('lab',comb),f'ms_{layer}_nominal':('lab05',comb),f'ms_{layer}_cswide':('lab_cw',comb),p+'opposite_sig':('opp_sig','any'),p+'entries':('study',lambda s:';'.join(sorted(s)))})
    best=I.sort_values(['a','power'],ascending=[True,False]).drop_duplicates('a').set_index('a')[['variant_id','b_s','standard_error','p_value','power','z_exp']]
    best.columns=[p+'index_variant',p+'beta_aligned',p+'se',p+'p',p+'power',p+'zexp']
    E=agg.join(best).reset_index()
    info=dict(offset=OFF[layer],n_files=len(FILES[layer]),K_pairs=int(K),alpha=float(alpha),n_cs_assessed=int(len(E)),labels=E[f'ms_{layer}_primary'].value_counts(dropna=False).to_dict(),opp_sig=int(E[p+'opposite_sig'].sum()))
    return E,I,info
res={}
Egs,Igs,res['info_gs']=labels('gs'); Ehb,Ihb,res['info_hb']=labels('hb')
Egsm,_,res['info_gs_multi']=labels('gs',multi=True)
Egs['gs_multi_gene_entry']=False
Egs.to_parquet(A4+'ms_gs_labels.parquet',index=False); Ehb.to_parquet(A4+'ms_hb_labels.parquet',index=False); Igs.to_parquet(A4+'ms_gs_pairs.parquet',index=False); Ihb.to_parquet(A4+'ms_hb_pairs.parquet',index=False)
M=M.merge(Egs,on='a',how='left').merge(Ehb,on='a',how='left').merge(Egsm[['a','ms_gs_primary']].rename(columns={'ms_gs_primary':'ms_gs_primary_multi'}),on='a',how='left')
hbgenes=set(HBmap.gene_id); M['hb_evaluable']=M.ga.isin(hbgenes)
# ---- verdict helper
def verdict(eff,thr,direction,p_holm,stop):
    if stop: return '不确定'
    met = (eff>=thr) if direction=='ge' else (eff<=thr)
    right= (eff>1) if direction=='ge' else (eff<1)
    if met and p_holm<0.05: return '成立'
    if (met and right) or (p_holm<0.05 and right): return '部分成立'
    return '不成立'
# ---- A1 frozen model
Fall=H4mod.features(M); tr=M[M.ms_ext_primary.notna()&(M.train==1)]; ytr=tr.ms_ext_primary.astype(int).values
A1=json.load(open(R+'h3h4_ext.json')); a1coef=A1['H4_ext']['coef_full_std']
models={}
for nm,cols in [('full',H4mod.FULL),('B0',H4mod.B0),('B1',H4mod.B1)]:
    _,lr,(med,mu,sd)=H4mod.fit_pred(Fall.loc[tr.index],ytr,Fall.loc[tr.index[:2]],cols); models[nm]=(cols,lr,med,mu,sd)
coef_diff=max(abs(float(v)-a1coef[k]) for k,v in zip(H4mod.FULL,models['full'][1].coef_[0]))
res['A1_model_reproduction']=dict(n_train=int(len(tr)),max_abs_coef_diff_vs_A1=coef_diff,ok=bool(coef_diff<1e-6))
def score(nm,idx):
    cols,lr,med,mu,sd=models[nm]; X=Fall.loc[idx,cols].fillna(med); return lr.predict_proba((X-mu)/sd)[:,1]
def h4_layer(lab,subset=None):
    te=M[M[lab].notna()&(M.train==0)] if subset is None else M[M[lab].notna()&subset]
    y=te[lab].astype(int).values; o=dict(n_test=int(len(te)),pos_test=int(y.sum()))
    o['stop_rule_triggered']=bool(len(te)<40 or min(y.sum(),len(y)-y.sum())<10)
    if len(te)<5 or y.min()==y.max(): return o
    pf,p0,p1=score('full',te.index),score('B0',te.index),score('B1',te.index)
    o['delong_full_vs_B0']=H4mod.delong(y,pf,p0); o['delong_full_vs_B1']=H4mod.delong(y,pf,p1)
    o['auc_full']=o['delong_full_vs_B0']['auc1']; o['auc_B0']=o['delong_full_vs_B0']['auc2']; o['auc_B1']=o['delong_full_vs_B1']['auc2']
    return o
def layer_block(layer):
    lab=f'ms_{layer}_primary'; out={}
    out['H3'],T=H3mod.h3_block(M,lab); n_det=out['H3']['n_evaluable']; out['layer_stop_lt100']=bool(n_det<100)
    out['H4']=h4_layer(lab); out['H4_allchr']=h4_layer(lab,subset=M.index.isin(M.index))
    p=dict(H3i=out['H3']['H3i_logit']['A_eqtl']['P'] if 'A_eqtl' in out['H3']['H3i_logit'] else 1.0,H3ii=out['H3']['H3ii_CNP_vs_CA']['P'],H4=out['H4'].get('delong_full_vs_B0',{}).get('P_one_sided',1.0))
    out['holm']=holm({f'H3i-{layer}':p['H3i'],f'H3ii-{layer}':p['H3ii'],f'H4-{layer}':p['H4']})
    return out,T
for layer in ['gs','hb']:
    res[layer],T=layer_block(layer); T.to_csv(A4+f'h3_{layer}_table.csv',index=False)
# ---- secondary
lp=M.log10p.abs(); zz=norm.isf(np.minimum(10**(-lp.clip(upper=300))/2,0.5)); zz=np.where(lp>=300,40,zz)
M['flag_beta_inconsistent']=((M.beta.abs()/(zz/np.sqrt(2*M.maf*(1-M.maf)*34000)))>5).astype(int)
sec={}
sec_flag_n=int(M.flag_beta_inconsistent.sum())
for layer in ['gs','hb']:
    lab=f'ms_{layer}_primary'; T=M[M[lab].notna()]
    sec[f'{layer}_PAV_fisher']=fisher(T.assign(y=T[lab].astype(int)),'PAV','y')
    sec[f'{layer}_nominal_H3'],_=H3mod.h3_block(M,f'ms_{layer}_nominal',allow_stop=False)
    sec[f'{layer}_cswide_H3'],_=H3mod.h3_block(M,f'ms_{layer}_cswide',allow_stop=False)
    sec[f'{layer}_excl_beta_flag_H3'],_=H3mod.h3_block(M[M.flag_beta_inconsistent==0] if 'flag_beta_inconsistent' in M else M,lab,allow_stop=False)
    M['tier']=H4mod.tiers(M); d=M[M[lab].notna()&(M.tier>0)]
    sec[f'{layer}_tiers']=dict(rates=d.groupby('tier')[lab].agg(['sum','size','mean']).reset_index().to_dict('records'),trend=H4mod.cochran_armitage(d.tier,d[lab]))
sec['gs_multi_gene_H3'],_=H3mod.h3_block(M,'ms_gs_primary_multi',allow_stop=False)
new=M[M.seer_ext_studies.isna()]; sec['gs_new_genes_H3'],_=H3mod.h3_block(new,'ms_gs_primary',allow_stop=False); sec['hb_new_genes_H3'],_=H3mod.h3_block(new,'ms_hb_primary',allow_stop=False)
def kap(a,b):
    d=M[M[a].notna()&M[b].notna()]
    if len(d)<2 or d[a].nunique()<2 and d[b].nunique()<2: return dict(n=int(len(d)))
    return dict(n=int(len(d)),agree=float((d[a]==d[b]).mean()),kappa=float(cohen_kappa_score(d[a],d[b])),table=pd.crosstab(d[a],d[b]).values.tolist())
sec['kappa']={f'{a}_vs_{b}':kap(a,b) for a,b in [('ms_gs_primary','ms_ext_primary'),('ms_gs_primary','ms_primary'),('ms_gs_primary','ms_hb_primary'),('ms_hb_primary','ms_ext_primary'),('ms_hb_primary','ms_primary')]}
def combl(r):
    v=[r.ms_gs_primary,r.ms_hb_primary]; v=[x for x in v if pd.notna(x)]
    return 1.0 if 1.0 in v else (0.0 if 0.0 in v else np.nan)
M['ms_independent_combined']=M.apply(combl,axis=1)
sec['combined_H3'],_=H3mod.h3_block(M,'ms_independent_combined',allow_stop=False)
# S-PAVdir
pav=con.sql(f"select * from '{B}pqtl_cs_pav_alltx_vep116.parquet' limit 1").df(); sec['pav_table_cols']=list(pav.columns)
PV=con.sql(f"select * from '{B}pqtl_cs_pav_alltx_vep116.parquet'").df()
PV=PV[PV.cs95_pav_cognate_alltx_vep116==1][['studyLocusId','cs95_pav_variants_alltx_vep116']].rename(columns={'studyLocusId':'a','cs95_pav_variants_alltx_vep116':'pav_variant'})
import re as _re
PV['pav_variant']=PV.pav_variant.apply(lambda v: _re.findall(r'(?:\d+|X)_\d+_[ACGT]+_[ACGT]+',str(v)))
PV=PV.explode('pav_variant')
sec['pav_variant_rows']=int(len(PV))
PV=PV.merge(L0[['a','variant_id','pp','b_ukb_locus']].rename(columns={'variant_id':'pav_variant'}),on=['a','pav_variant']).sort_values(['a','pp'],ascending=[True,False]).drop_duplicates('a')
M=M.merge(PV[['a','pav_variant','b_ukb_locus']].rename(columns={'b_ukb_locus':'pav_ukb_beta'}),on='a',how='left')
M['ms_gs_pavdir']=np.where((M.PAV==1)&(M.ms_gs_primary==1)&(M.pav_ukb_beta<0),np.nan,M.ms_gs_primary)
sec['S_PAVdir']=dict(n_relabelled=int(((M.PAV==1)&(M.ms_gs_primary==1)&(M.pav_ukb_beta<0)).sum()),n_pav_supported=int(((M.PAV==1)&(M.ms_gs_primary==1)).sum()),n_pav_supported_no_pavbeta=int(((M.PAV==1)&(M.ms_gs_primary==1)&M.pav_ukb_beta.isna()).sum()))
sec['S_PAVdir_H3'],_=H3mod.h3_block(M,'ms_gs_pavdir',allow_stop=False)
T=M[M.ms_gs_pavdir.notna()]; sec['S_PAVdir_PAV_fisher']=fisher(T.assign(y=T.ms_gs_pavdir.astype(int)),'PAV','y')
# HOLBAEK peptide tier sensitivity (closest analogue: count support only when ST4 cis pQTL in CS is tier1)
st4=pd.read_csv(A4+'niu_st4_tiers.tsv',sep='\t'); st4=st4[st4.cistrans=='cis']
rs=con.sql(f"select variant_id, unnest(rsIds) rsid from '{B}variant_core.parquet'").df(); rs=rs[rs.rsid.isin(set(st4.rsID))]
st4=st4.merge(rs,left_on='rsID',right_on='rsid')
lt=L0[['a','variant_id']].merge(st4[['variant_id','tier']],on='variant_id').groupby('a').tier.apply(lambda s:'tier1' if (s=='tier1').any() else sorted(s)[0])
M['hb_st4_tier']=M.a.map(lt)
M['ms_hb_tier1']=np.where((M.ms_hb_primary==1)&(M.hb_st4_tier!='tier1'),np.nan,M.ms_hb_primary)  # strict: support only if ST4 tier1 cis pQTL in CS
M['ms_hb_tierflag']=np.where((M.ms_hb_primary==1)&M.hb_st4_tier.isin(['tier2','tier3','tier4']),np.nan,M.ms_hb_primary)  # relabel only single-peptide-tier CS
sec['HB_tier']=dict(cs_with_st4=int(M.hb_st4_tier.notna().sum()),tiers=M.hb_st4_tier.value_counts().to_dict(),relabelled_strict=int(((M.ms_hb_primary==1)&(M.hb_st4_tier!='tier1')).sum()),relabelled_flag=int(((M.ms_hb_primary==1)&M.hb_st4_tier.isin(['tier2','tier3','tier4'])).sum()))
sec['HB_tierflag_H3'],_=H3mod.h3_block(M,'ms_hb_tierflag',allow_stop=False)
T=M[M.ms_hb_tierflag.notna()]; sec['HB_tierflag_PAV_fisher']=fisher(T.assign(y=T.ms_hb_tierflag.astype(int)),'PAV','y')
sec['HB_tier1_H3'],_=H3mod.h3_block(M,'ms_hb_tier1',allow_stop=False)
T=M[M.ms_hb_tier1.notna()]; sec['HB_tier1_PAV_fisher']=fisher(T.assign(y=T.ms_hb_tier1.astype(int)),'PAV','y')
sec['n_flag_beta_inconsistent']=sec_flag_n
res['secondary']=sec
# verdicts (GS confirmatory)
g=res['gs']; stop=g['layer_stop_lt100']
h3i=g['H3']['H3i_logit'].get('A_eqtl',{}); h3ii=g['H3']['H3ii_CNP_vs_CA']; h4=g['H4']
res['verdicts']={'H3i-A4':verdict(h3i.get('OR',np.nan),1.5,'ge',g['holm']['H3i-gs'],stop),
 'H3ii-A4':verdict(h3ii.get('OR',np.nan),0.5,'le',g['holm']['H3ii-gs'],stop or g['H3']['H3ii_stop_rule_triggered']),
 'H4-A4':verdict(1+h4.get('delong_full_vs_B0',{}).get('dauc',np.nan),1.05,'ge',g['holm']['H4-gs'],stop or h4['stop_rule_triggered'])}
h=res['hb']; stoph=h['layer_stop_lt100']
res['verdicts_hb']={'H3i-HB':verdict(h['H3']['H3i_logit'].get('A_eqtl',{}).get('OR',np.nan),1.5,'ge',h['holm']['H3i-hb'],stoph),
 'H3ii-HB':verdict(h['H3']['H3ii_CNP_vs_CA'].get('OR',np.nan),0.5,'le',h['holm']['H3ii-hb'],stoph or h['H3']['H3ii_stop_rule_triggered']),
 'H4-HB':verdict(1+h['H4'].get('delong_full_vs_B0',{}).get('dauc',np.nan),1.05,'ge',h['holm']['H4-hb'],stoph or h['H4']['stop_rule_triggered'])}
keep=['a','ms_gs_primary','ms_gs_nominal','ms_gs_cswide','gs_opposite_sig','gs_entries','gs_index_variant','gs_beta_aligned','gs_se','gs_p','gs_power','gs_zexp','ms_gs_primary_multi',
 'ms_hb_primary','ms_hb_nominal','ms_hb_cswide','hb_opposite_sig','hb_entries','hb_index_variant','hb_beta_aligned','hb_se','hb_p','hb_power','hb_zexp','hb_evaluable','hb_st4_tier','ms_independent_combined','ms_gs_pavdir','ms_hb_tier1','ms_hb_tierflag','pav_variant','pav_ukb_beta']
M[keep].to_parquet(A4+'a4_cs_fields.parquet',index=False)
dump(res,'a4.json'); print(json.dumps({k:res[k] for k in ['info_gs','info_hb','A1_model_reproduction','verdicts','verdicts_hb']},default=str))
