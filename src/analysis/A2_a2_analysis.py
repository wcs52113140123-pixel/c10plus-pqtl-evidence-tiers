"""Addendum A2 analysis: H2a/H2b reproduction in ARIC (EA, AA) and Fenland (SomaScan v4). Prereg A2 SHA256 08df92f0..."""
import sys, json, numpy as np, pandas as pd, duckdb
sys.path.insert(0,'.')
from c10common import *
from scipy.stats import pearsonr, fisher_exact, norm
import statsmodels.formula.api as smf
B="../data_v1/db/parquet/"
M=derive(load())
L=pd.read_parquet(W+'ukb_cis_locus_grch37.parquet'); L['chrom']=L.chrom.astype(str); L['cp']=list(zip(L.chrom,L.pos_b38.astype(int)))
p95=L[L.is95].groupby('a').cp.apply(set).to_dict(); p99=L.groupby('a').cp.apply(set).to_dict()
vinfo={(a,cp):(v,r,al,b) for a,cp,v,r,al,b in zip(L.a,L.cp,L.variant_id,L.ref,L.alt,L.beta)}
lead=M.set_index('a')[['lead_variant_id','beta']].to_dict('index')
def ukb_beta(a,v,b):
    return lead[a]['beta'] if lead[a]['lead_variant_id']==v else b
def build(S,tag,strict=True):
    """S: gene_id, chr, pos, src, pip, allele, beta, p. Returns per-cs columns."""
    S=S.dropna(subset=['gene_id','pos']).copy(); S['chr']=S.chr.astype(str).str.replace(r'\.0$','',regex=True).replace({'23':'X'}); S['cp']=list(zip(S.chr,S.pos.astype(int)))
    Sg={g:x for g,x in S.groupby('gene_id')}
    rows=[]
    for a,g in zip(M.a,M.ga):
        x=Sg.get(g); r=dict(a=a)
        r[f'{tag}_has_cis']=int(x is not None)
        if x is None: rows.append(r); continue
        m95=x[[cp in p95.get(a,set()) for cp in x.cp]]
        r[f'rep_{tag}']=int(len(m95)>0)
        r[f'rep_{tag}_99']=int(any(cp in p99.get(a,set()) for cp in x.cp))
        if strict: r[f'rep_{tag}_strict']=int(((m95.src=='finemap')&(m95.pip>=0.5)).any())
        d=m95[m95.beta.notna()&m95.allele.notna()].sort_values('p')
        for q in d.itertuples():
            v,ref,alt,bl=vinfo[(a,q.cp)]
            sg=1 if q.allele==alt else (-1 if q.allele==ref else 0)
            if sg: r[f'{tag}_b_ukb']=ukb_beta(a,v,bl); r[f'{tag}_b_rep']=sg*q.beta; break
        rows.append(r)
    return pd.DataFrame(rows)
def block(D,rep):
    o={'n_cs':int(len(D)),'n_genes':int(D.ga.nunique()),'n_rep':int(D[rep].sum()),'rep_rate':float(D[rep].mean())}
    o['stop_rule_triggered']=bool(D[rep].sum()<30)
    cov=['mlp','logn95','PAV']
    o['main'],_=glm(D,rep,['A_gtex',*cov]); o['negctrl_wrong'],_=glm(D,rep,['W_gtex',*cov])
    j,mc=glm(D,rep,['A_gtex','W_gtex',*cov],terms=['A_gtex','W_gtex']); o['joint']=j; o['joint_diff']=diff_test(mc,'A_gtex','W_gtex')
    o['nongtex'],_=glm(D,rep,['A_nongtex',*cov])
    S_=D[D.mlog10p>30]; o['power_matched_subset'],_=glm(S_,rep,['A_gtex','mlp','logn95','PAV'])
    o['PAV_adj'],_=glm(D,rep,['PAV','mlp','logn95','A_gtex'])
    o['four_group']=D.groupby(['A_gtex','PAV'])[rep].agg(['mean','size']).reset_index().to_dict('records')
    D=D.assign(_s=1)
    Mt=pav_match(D,rep,'_s')
    o['pav_matched_n_sets']=int(Mt.set_id.nunique()) if Mt is not None else 0; o['pav_n_treated_total']=int(D.PAV.sum())
    o['H2b_stop_rule_triggered']=bool(o['pav_matched_n_sets']<50)
    if Mt is not None:
        o['pav_matched'],_=glm(Mt,rep,['PAV'],cluster='set_id'); o['pav_matched_rates']={'PAV1':float(Mt[Mt.PAV==1][rep].mean()),'PAV0':float(Mt[Mt.PAV==0][rep].mean())}
        o['pav_smd_before']=smd(D.dropna(subset=['mlp','absb','maf','logn95']),'PAV'); o['pav_smd_after']=smd(Mt,'PAV')
    d=D.dropna(subset=['mlp','absb','maf','logn95'])
    try:
        m=smf.glm(f"{rep} ~ PAV + cr(mlp, df=4) + cr(absb, df=4) + cr(maf, df=4) + cr(logn95, df=4)",d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(d.ga)[0]})
        o['pav_spline_adj']=dict(OR=float(np.exp(m.params['PAV'])),P=float(m.pvalues['PAV']),n=int(len(d)))
    except Exception as e: o['pav_spline_adj']={'error':str(e)}
    return o
def effcorr(D,tag):
    if f'{tag}_b_ukb' not in D: return {'n_pairs':0}
    E=D.dropna(subset=[f'{tag}_b_ukb',f'{tag}_b_rep']).copy(); E=E[(E[f'{tag}_b_ukb']!=0)&(E[f'{tag}_b_rep']!=0)]
    E['opp']=(np.sign(E[f'{tag}_b_ukb'])!=np.sign(E[f'{tag}_b_rep'])).astype(int); o={'n_pairs':int(len(E))}
    for k,gp in E.groupby('PAV'):
        o[f'PAV{k}']=dict(n=int(len(gp)),r=float(pearsonr(gp[f'{tag}_b_ukb'],gp[f'{tag}_b_rep'])[0]) if len(gp)>3 else None,opp_frac=float(gp.opp.mean()),n_opp=int(gp.opp.sum()))
    if 'PAV0' in o and 'PAV1' in o and o['PAV1']['r'] is not None:
        o['fisher_z_PAV1_lt_PAV0']=fisher_z_diff(o['PAV1']['r'],o['PAV1']['n'],o['PAV0']['r'],o['PAV0']['n'])
        t=[[o['PAV1']['n_opp'],o['PAV1']['n']-o['PAV1']['n_opp']],[o['PAV0']['n_opp'],o['PAV0']['n']-o['PAV0']['n_opp']]]
        OR,p=fisher_exact(t,alternative='greater'); o['opp_fisher_PAV1_gt']=dict(table=t,OR=float(OR),P_one_sided=float(p))
    return o
res={}; cols=[]
for anc in ['EA','AA']:
    S=pd.read_csv(W+f'a2_aric_signals_{anc}.tsv',sep='\t',dtype={'chr':str}); tag=f'aric_{anc.lower()}'
    C=build(S,tag); cols.append(C); M=M.merge(C,on='a',how='left')
# Fenland
A=pd.read_csv(W+'a2_fenland_assoc_mapped.tsv',sep='\t',dtype={'chr':str})
A=A[A.cis==True].copy()
miss=A.pos.isna()
if miss.any():
    con=duckdb.connect(); rm=con.sql(f"select r rsid, chrom, pos_b38 from (select unnest(rsIds) r, chrom, pos_b38 from '{B}variant_core.parquet') where r in (select unnest(?::varchar[]))",params=[list(A.loc[miss,'rsid'].dropna().astype(str))]).df().drop_duplicates('rsid')
    A=A.merge(rm,on='rsid',how='left'); A.loc[miss,'chr']=A.loc[miss,'chrom_y'] if 'chrom_y' in A else A.loc[miss,'chrom']; A.loc[miss,'pos']=A.loc[miss,'pos_b38']
pv=pd.to_numeric(A.p_value,errors='coerce')
bnum=pd.to_numeric(A.beta.astype(str).str.extract(r'([0-9.]+)')[0],errors='coerce'); bsign=np.where(A.beta.astype(str).str.contains('decrease'),-1.0,np.where(A.beta.astype(str).str.contains('increase'),1.0,np.nan))
SF=pd.DataFrame({'gene_id':A.gene_id,'chr':A.chr.astype(str),'pos':A.pos,'src':'lead','pip':np.nan,'allele':A.effect_allele,'beta':bnum*bsign,'p':pv})
res['fenland_beta_parsed']=int(np.isfinite(bnum*bsign).sum())
C=build(SF,'fenland',strict=False); M=M.merge(C,on='a',how='left')
res['fenland_cis_assoc_used']=int(len(SF)); res['fenland_pos_missing_filled_by_rsid']=int(miss.sum())
for tag in ['aric_ea','aric_aa','fenland']:
    D=M[M[f'{tag}_has_cis']==1].copy()
    res[tag]=block(D,f'rep_{tag}')
    res[tag]['effect_corr']=effcorr(D,tag)
    res[tag]['sens_rep99'],_=glm(D,f'rep_{tag}_99',['A_gtex','mlp','logn95','PAV'])
    if f'rep_{tag}_strict' in D:
        res[tag]['sens_strict'],_=glm(D,f'rep_{tag}_strict',['A_gtex','mlp','logn95','PAV'])
        Mt=pav_match(D.assign(_s=1),f'rep_{tag}_strict','_s'); res[tag]['sens_strict_pav_matched'],_=glm(Mt,f'rep_{tag}_strict',['PAV'],cluster='set_id')
    D8=D[D.v4&~D.v3]
    res[tag]['sens_H2genes_main'],_=glm(D8,f'rep_{tag}',['A_gtex','mlp','logn95','PAV'])
    Mt=pav_match(D8.assign(_s=1),f'rep_{tag}','_s'); res[tag]['sens_H2genes_pav_matched'],_=glm(Mt,f'rep_{tag}',['PAV'],cluster='set_id')
    res[tag]['overlap_genes_with_deCODE_v4']=float(D.v4.mean())
# confirmatory family Holm
P={'H2a-ARIC':res['aric_ea']['main']['A_gtex']['P'],'H2b-ARIC':res['aric_ea']['pav_matched']['PAV']['P'],'H2a-Fenland':res['fenland']['main']['A_gtex']['P'],'H2b-Fenland':res['fenland']['pav_matched']['PAV']['P']}
res['holm_A2']=holm(P); res['raw_P_A2']=P
# random-effects meta with deCODE (H2)
h12=json.load(open(R+'h1_h2.json'))
def dl(ests):
    y=np.array([np.log(e['OR']) for e in ests]); se=np.array([(np.log(e['hi'])-np.log(e['lo']))/3.92 for e in ests]); w=1/se**2
    yf=(w*y).sum()/w.sum(); Q=(w*(y-yf)**2).sum(); k=len(y); C_=w.sum()-(w**2).sum()/w.sum(); tau2=max(0,(Q-(k-1))/C_)
    ws=1/(se**2+tau2); yr=(ws*y).sum()/ws.sum(); ser=ws.sum()**-0.5
    return dict(OR=float(np.exp(yr)),lo=float(np.exp(yr-1.96*ser)),hi=float(np.exp(yr+1.96*ser)),P=float(2*norm.sf(abs(yr/ser))),tau2=float(tau2),I2=float(max(0,(Q-(k-1))/Q)) if Q>0 else 0.0,Q=float(Q),k=k)
for nm,key,term in [('H2a','main','A_gtex'),('H2b','pav_matched','PAV')]:
    e_dec=h12['H2'][key][term]; e_ar=res['aric_ea'][key][term]; e_fe=res['fenland'][key][term]
    res[f'meta_{nm}_3cohort']=dl([e_dec,e_ar,e_fe]); res[f'meta_{nm}_ARIC_Fenland']=dl([e_ar,e_fe])
keep=['a']+[c for c in M.columns if c.startswith(('aric_','rep_aric','fenland_','rep_fenland'))]
M[keep].to_parquet(W+'a2_cs_fields.parquet',index=False)
dump(res,'a2.json')
for tag in ['aric_ea','aric_aa','fenland']:
    x=res[tag]; print(tag,x['n_cs'],x['n_genes'],x['n_rep'],round(x['rep_rate'],3),'H2a',{k:round(v,4) if isinstance(v,float) else v for k,v in x['main']['A_gtex'].items() if k in ('OR','lo','hi','P')},'diffP',round(x['joint_diff']['P_one_sided'],6),'H2b',{k:round(v,4) if isinstance(v,float) else v for k,v in x.get('pav_matched',{}).get('PAV',{}).items() if k in ('OR','lo','hi','P','n')},'sets',x['pav_matched_n_sets'])
print(res['holm_A2']); print(res['meta_H2a_3cohort'],res['meta_H2b_3cohort'])
