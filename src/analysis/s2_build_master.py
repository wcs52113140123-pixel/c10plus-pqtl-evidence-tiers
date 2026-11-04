"""Stage 2: build the per-credible-set master table (exposures, covariates, replication, truth layers).
Implements definitions frozen in C10plus_preregistration_v2_frozen.md (SHA256 126f515f...)."""
import duckdb, numpy as np, pandas as pd, re, time, json
from scipy.special import logsumexp
t0=time.time()
B="../data_v1/db/parquet/"
W="./work/"
con=duckdb.connect(); con.sql("SET threads=96; SET memory_limit='150GB'")
log={}
# ---------------- coloc from sufficient statistics ----------------
P=pd.read_parquet(W+"pairs_stats.parquet")
def h4(P,p1=1e-4,p2=1e-4,p12=1e-5,mode='union'):
    k=P.k.values.astype(float); na=P.na.values.astype(float); nb=P.nb.values.astype(float)
    Sa=P.Sa.values; Sb=P.Sb.values; L12=P.Lsh12.values; L1=P.Lsh1.values; L2=P.Lsh2.values
    if mode=='union':
        with np.errstate(divide='ignore',invalid='ignore'):
            s1=np.logaddexp(Sa,np.log(np.maximum(nb-k,0)))
            s2=np.logaddexp(Sb,np.log(np.maximum(na-k,0)))
            d1=np.minimum(L1-Sa,0); d2=np.minimum(L2-Sb,0)
            aonly=np.where(d1>-1e-9,-np.inf,Sa+np.log1p(-np.exp(d1)))
            bonly=np.where(d2>-1e-9,-np.inf,Sb+np.log1p(-np.exp(d2)))
            s12=np.logaddexp(np.logaddexp(L12,aonly),bonly)
    else:
        s1,s2,s12=L1,L2,L12
    r=np.minimum(np.exp(s12-(s1+s2)),1-1e-12)
    lH=np.vstack([np.zeros_like(s1),np.log(p1)+s1,np.log(p2)+s2,np.log(p1)+np.log(p2)+s1+s2+np.log1p(-r),np.log(p12)+s12])
    return np.exp(lH[4]-logsumexp(lH,axis=0))
P['h4']=h4(P); P['h4_p12_5e6']=h4(P,p12=5e-6); P['h4_p12_1e6']=h4(P,p12=1e-6); P['h4_inter']=h4(P,mode='inter')
st=P.studyType; pj=P.projectId; cog=P.gb==P.ga
P['cat']=np.select([pj.eq('Sun_2018'),(st.eq('eqtl')&pj.eq('GTEx-v10')&cog),(st.isin(['eqtl','sceqtl'])&~pj.eq('GTEx-v10')&cog),(st.isin(['sqtl','tuqtl'])&cog),(st.eq('eqtl')&pj.eq('GTEx-v10')&~cog),P.ukbtrans.fillna(False)],
                   ['sun','gtex','nongtex','iso','wrong','ukbtrans'],'other')
log['pairs_by_cat']=P['cat'].value_counts().to_dict()
U=pd.read_parquet(W+"ukb_cis_cs.parquet")
M=U.copy()
for cat in ['sun','gtex','nongtex','iso','wrong']:
    g=P[P['cat']==cat].groupby('a')
    agg=g[['h4','h4_p12_5e6','h4_p12_1e6','h4_inter']].max().add_prefix(cat+'_')
    M=M.merge(agg,left_on='a',right_index=True,how='left')
for c_ in [c for c in M.columns if '_h4' in c]: M[c_]=M[c_].fillna(0.0)
# number of anchoring eQTL datasets (distinct studyId with h4>=0.8)
an=P[P['cat'].isin(['gtex','nongtex'])&(P.h4>=0.8)].groupby('a').studyId.nunique().rename('n_anchor_datasets')
M=M.merge(an,left_on='a',right_index=True,how='left'); M['n_anchor_datasets']=M.n_anchor_datasets.fillna(0)
# Sun direction at top shared variant of best Sun pair
sb=P[P['cat']=='sun'].sort_values('h4',ascending=False).drop_duplicates('a')[['a','b1top','b2top','vtop']].rename(columns={'b1top':'sun_b_ukb','b2top':'sun_b_sun','vtop':'sun_vtop'})
M=M.merge(sb,on='a',how='left')
# hotspot: distinct other genes with UKB trans CS sharing a variant
hs=P[(P['cat']=='ukbtrans')&(P.gb!=P.ga)].groupby('a').gb.nunique().rename('n_trans_genes_shared')
M=M.merge(hs,left_on='a',right_index=True,how='left'); M['n_trans_genes_shared']=M.n_trans_genes_shared.fillna(0)
print('coloc done',time.time()-t0,flush=True)
# ---------------- UKB attributes ----------------
A=con.sql(f"""select c.studyLocusId a, c.studyId, c.lead_variant_id, c.chrom, c.pos_b38, c.beta, c.eaf, c.log10p, c.purityMinR2, c.n_locus,
   g.tss, g.symbol, s.assay_id oid from '{B}qtl_cs.parquet' c left join '{B}gene.parquet' g using(gene_id) left join '{B}pqtl_study.parquet' s using(studyId)
   where c.projectId='UKB_PPP_EUR' and not c.isTransQtl""").df()
L=con.sql(f"""select l.studyLocusId a, l.variant_id v, l.pp, l.beta, l.is95, vc.allele_key, vc.ref, vc.alt, vc.chrom vchrom, vc.pos_b38 vpos,
   list_filter(vc.alleleFrequencies, x -> x.populationName='nfe_adj')[1].alleleFrequency nfe_af
   from read_parquet('{B}qtl_cs_locus/*/*.parquet', hive_partitioning=true) l join '{B}variant_core.parquet' vc using(variant_id)
   where l.studyType='pqtl' and l.studyLocusId in (select a from A)""").df()
log['ukb_locus_rows']=len(L); log['ukb_locus_rows_missing_vc']=int(con.sql(f"select count(*) from read_parquet('{B}qtl_cs_locus/*/*.parquet', hive_partitioning=true) l where l.studyType='pqtl' and l.studyLocusId in (select a from A)").fetchone()[0])-len(L)
n95=L.groupby('a').is95.sum().rename('n95')
A=A.merge(n95,left_on='a',right_index=True,how='left')
M=M.merge(A,on='a',how='left')
lf=L[['a','v','nfe_af']].rename(columns={'v':'lead_variant_id','nfe_af':'lead_nfe_af'})
M=M.merge(lf,on=['a','lead_variant_id'],how='left')
pav=con.sql(f"""select p.studyLocusId a, a2.cs95_pav_cognate_alltx_vep116 pav, a2.pp_sum_pav_cognate_alltx_vep116 pav_pp, p.cs95_pav_cognate_canonical pav_canon, lp.lead_is_pav_cognate pav_lead
   from '{B}pqtl_cs_pav.parquet' p join '{B}pqtl_cs_pav_alltx_vep116.parquet' a2 using(studyLocusId) left join '{B}pqtl_cs_leadpav.parquet' lp using(studyLocusId)""").df()
M=M.merge(pav,on='a',how='left')
cov=con.sql(f"select gene_id ga, assayed_SomaScan_v3_Sun2018 v3, assayed_SomaScan_v4_deCODE v4, assayed_MS_Seer ms_assayed from '{B}coverage_ukbppp_cis_gene_matrix.parquet'").df()
M=M.merge(cov,on='ga',how='left')
sungenes=set(con.sql(f"select distinct gene_id from '{B}qtl_cs.parquet' where projectId='Sun_2018' and not isTransQtl").df().gene_id)
M['sun_has_cis']=M.ga.isin(sungenes).astype(int)
# ---------------- deCODE matching ----------------
D=con.sql(f"""select layer, assay_id, gene_id, chrom, pos_b38, effect_allele, other_allele, allele_key, beta, mlog10p, signal_id, cs_members_src
   from '{B}ext_pqtl_layers.parquet' where layer in ('decode_somascan_ST19','decode_somascan_ST33_cs') and (cis_trans_src='cis' or is_cis_1mb) and gene_id is not null""").df()
M['decode_has_cis']=M.ga.isin(set(D.gene_id)).astype(int)
rs=set(); 
for s in D.cs_members_src.dropna(): rs.update(t for t in s.split(',') if t.startswith('rs'))
rsmap=con.sql(f"select r rsid, chrom, pos_b38 from (select unnest(rsIds) r, chrom, pos_b38 from '{B}variant_core.parquet') where r in (select unnest(?::varchar[]))",params=[list(rs)]).df()
rsd=rsmap.groupby('rsid').apply(lambda x: list(zip(x.chrom.astype(str),x.pos_b38.astype(int)))).to_dict()
def members(r):
    out={(str(r.chrom),int(r.pos_b38))} if pd.notna(r.pos_b38) else set()
    if isinstance(r.cs_members_src,str):
        for t in r.cs_members_src.split(','):
            m=re.match(r'chr(\w+):(\d+)$',t)
            if m: out.add((m.group(1),int(m.group(2))))
            elif t in rsd: out.update(rsd[t])
    return out
D['mem']=[members(r) for r in D.itertuples()]
D['leadpos']=[(str(r.chrom),int(r.pos_b38)) for r in D.itertuples()]
L['cp']=list(zip(L.vchrom.astype(str),L.vpos.astype(int)))
pos95=L[L.is95].groupby('a').cp.apply(set).to_dict(); pos99=L.groupby('a').cp.apply(set).to_dict()
ak99=L.groupby('a').apply(lambda x: dict(zip(x.allele_key,zip(x.beta,x.ref,x.alt,x.v)))).to_dict()
Dg={g:x for g,x in D.groupby('gene_id')}
rows=[]
for r in M[['a','ga']].itertuples():
    x=Dg.get(r.ga); p95=pos95.get(r.a,set()); p99=pos99.get(r.a,set())
    rec=dict(a=r.a,rep_decode=0,rep_decode_leadonly=0,rep_decode_99=0)
    if x is not None:
        m95=[i for i,mm in enumerate(x.mem) if mm & p95]
        rec['rep_decode']=int(len(m95)>0)
        rec['rep_decode_99']=int(any(mm & p99 for mm in x.mem))
        rec['rep_decode_leadonly']=int(any(lp in p99 for lp in x.leadpos))
        # direction: matched signals (95 rule) with lead allele in UKB locus; take max mlog10p
        cand=x.iloc[m95] if m95 else x.iloc[[i for i,mm in enumerate(x.mem) if mm & p99]]
        cand=cand.sort_values('mlog10p',ascending=False)
        akd=ak99.get(r.a,{})
        for q in cand.itertuples():
            if isinstance(q.allele_key,str) and q.allele_key in akd:
                bu,ref,alt,v=akd[q.allele_key]
                sgn=1 if q.effect_allele==alt else (-1 if q.effect_allele==ref else 0)
                if sgn!=0 and pd.notna(q.beta):
                    rec.update(dec_b_ukb=bu,dec_b_dec=sgn*q.beta,dec_v=v,dec_mlp=q.mlog10p); break
    rows.append(rec)
M=M.merge(pd.DataFrame(rows),on='a',how='left')
print('deCODE done',time.time()-t0,flush=True)
# ---------------- MS: MSPA (Suhre 2025 ST10/ST9) ----------------
gm=con.sql(f"select gene_id, unnest(uniprot_all) u from '{B}gene.parquet'").df()
u2g=gm.groupby('u').gene_id.apply(set).to_dict()
X=pd.ExcelFile(W+"suhre2025_MSPA.xlsx")
ms=[]
for sh,repcol in [('ST10','replicated (p<0.05/392)'),('ST9','replicated (p<0.05/319)')]:
    d=pd.read_excel(X,sheet_name=sh,header=2) if False else pd.read_excel(X,sheet_name=sh)
    if sh=='ST10': chrom=d['CHROM'].astype(str); pos=pd.to_numeric(d['GENPOS (hg38)'],errors='coerce'); upc='Target UniProt'
    else:
        mm=d['SNP (build 38)'].astype(str).str.extract(r'(?:chr)?(\w+?)[:_](\d+)'); chrom=mm[0]; pos=pd.to_numeric(mm[1],errors='coerce'); upc='UniProt'
    for i in range(len(d)):
        ups=re.split(r'[;,| ]+',str(d[upc].iloc[i])); genes=set().union(*[u2g.get(u.split('-')[0],set()) for u in ups]) if ups else set()
        pav_rs=str(d['MBRplus.PEPanno.PAVrsid'].iloc[i]).strip()
        ms.append(dict(sheet=sh,row=i,rank=d['Rank'].iloc[i],chrom=str(chrom.iloc[i]).replace('chr',''),pos=pos.iloc[i],genes=genes,
            replicated=bool(d[repcol].iloc[i]),powered=bool(d['power>99%'].iloc[i]),power=float(pd.to_numeric(d['power'].iloc[i],errors='coerce')),
            mspa_hi=bool(d['MSPA>0.8'].iloc[i]),mspa_lo=bool(d['MSPA<0.2'].iloc[i]),pep_pav=pav_rs not in ('.','nan','','None')))
MS=pd.DataFrame(ms); log['mspa_rows']=MS.sheet.value_counts().to_dict()
Lg=L.merge(M[['a','ga']],on='a')
key=Lg.assign(chrom=Lg.vchrom.astype(str),pos=Lg.vpos.astype(int))[['a','ga','chrom','pos','pp']]
mrows=[]
for r in MS.itertuples():
    if pd.isna(r.pos) or not r.genes: continue
    k_=key[(key.chrom==r.chrom)&(key.pos==int(r.pos))&(key.ga.isin(r.genes))]
    if len(k_): mrows.append(dict(r._asdict(),a=k_.sort_values('pp',ascending=False).a.iloc[0]))
MM=pd.DataFrame(mrows)
log['mspa_matched_rows']=MM.sheet.value_counts().to_dict() if len(MM) else {}
MM['pri']=(MM.sheet=='ST10').astype(int)
MM=MM.sort_values(['a','pri','power'],ascending=[True,False,False]).drop_duplicates('a')
MM['ms_primary']=np.where(MM.replicated,1.0,np.where(MM.powered,0.0,np.nan))
MM['ms_mspa']=np.where(MM.mspa_hi,1.0,np.where(MM.mspa_lo,0.0,np.nan))
M=M.merge(MM[['a','sheet','rank','ms_primary','ms_mspa','pep_pav','power','powered','replicated']].rename(columns={'sheet':'ms_sheet','rank':'ms_rank','pep_pav':'ms_pep_pav','power':'ms_power','powered':'ms_powered','replicated':'ms_replicated'}),on='a',how='left')
# ---------------- MS: Seer discovery layer ----------------
S=con.sql(f"""select layer, gene_id, allele_key, mlog10p, cis_trans_src, is_cis_1mb from '{B}ext_pqtl_layers.parquet' where layer like 'seer_ms%' and gene_id is not null""").df()
seer_cov=set(S[S.layer=='seer_ms_ST1'].gene_id)
sig=S[((S.layer=='seer_ms_ST2_primary')&((S.cis_trans_src=='True')|S.is_cis_1mb))|((S.layer=='seer_ms_ST3')&S.is_cis_1mb)|((S.layer=='seer_ms_ST1')&((S.cis_trans_src=='True')|S.is_cis_1mb)&(S.mlog10p>=7.3))]
sigk=sig.dropna(subset=['allele_key']).groupby('gene_id').allele_key.apply(set).to_dict()
aks=L.groupby('a').allele_key.apply(set).to_dict()
M['ms_seer']=[ (1.0 if (sigk.get(g,set()) & aks.get(a,set())) else 0.0) if g in seer_cov else np.nan for a,g in zip(M.a,M.ga)]
print('MS done',time.time()-t0,flush=True)
# ---------------- non-European ancestries ----------------
con.register('Lak',L[['a','v','allele_key','ref','alt','pp','beta','nfe_af']].merge(M[['a','ga','oid','lead_variant_id','beta']].rename(columns={'beta':'lead_beta'}),on='a'))
NE=con.sql(f"""select l.a, l.v, l.pp, l.beta eur_b_locus, l.lead_beta, l.lead_variant_id, l.ref, l.alt, l.nfe_af, n.ancestry, n.allele0, n.allele1, n.a1freq, n.info ninfo, n.beta nb, n.se nse
   from Lak l join read_parquet('{B}ukbppp_nonEUR_cis_sumstats/*/*.parquet', hive_partitioning=true) n on n.gene_id=l.ga and n.allele_key=l.allele_key and n.oid=l.oid""").df()
log['nonEUR_rows_joined']=len(NE)
NE['sgn']=np.where(NE.allele1==NE.alt,1,np.where((NE.allele1==NE.ref)&(NE.allele0==NE.alt),-1,0))
NE['maf']=np.minimum(NE.a1freq,1-NE.a1freq)
NE=NE[(NE.sgn!=0)&(NE.maf>=0.01)&((NE.ninfo>=0.8)|NE.ninfo.isna())&(NE.nse>0)]
NE['is_lead']=(NE.v==NE.lead_variant_id)
NE['eur_b']=np.where(NE.is_lead,NE.lead_beta,NE.eur_b_locus)
NE=NE.sort_values(['a','ancestry','is_lead','pp'],ascending=[True,True,False,False])
def meta(E,tag):
    E=E.copy(); E['b_al']=E.sgn*E.nb*np.sign(E.eur_b); E['w']=1/E.nse**2
    E['nfe_maf']=np.minimum(E.nfe_af,1-E.nfe_af); E['mafdiff']=(E.maf-E.nfe_maf).abs()
    g=E.groupby('a')
    out=pd.DataFrame({'xa_b':g.apply(lambda x:(x.w*x.b_al).sum()/x.w.sum()),'xa_se':g.w.sum()**-0.5,'xa_nanc':g.size(),'xa_maxmafdiff':g.mafdiff.max(),'xa_eur_absb':g.eur_b.apply(lambda s: s.abs().mean())})
    from scipy.stats import norm
    out['xa_p']=2*norm.sf(np.abs(out.xa_b/out.xa_se)); out['xa_rep']=((out.xa_b>0)&(out.xa_p<0.05)).astype(int); out['xa_samedir']=(out.xa_b>0).astype(int)
    return out.add_suffix(tag)
E1=NE.drop_duplicates(['a','ancestry'])
XA=meta(E1,'')
XL=meta(NE[NE.is_lead].drop_duplicates(['a','ancestry']),'_leadonly')
# per-ancestry aligned estimates for AFR/CSA and slope analysis
E1=E1.assign(b_h=E1.sgn*E1.nb)
for anc,tag in [('African','afr'),('Central_South Asian','csa')]:
    e=E1[E1.ancestry==anc].set_index('a')
    XA['xa_rep_'+tag]=((np.sign(e.b_h)==np.sign(e.eur_b))&(2*__import__('scipy.stats',fromlist=['norm']).norm.sf(np.abs(e.b_h/e.nse))<0.05)).astype(float).reindex(XA.index)
E1[['a','ancestry','v','is_lead','eur_b','b_h','nse','maf','nfe_af']].to_parquet(W+"nonEUR_index_variants.parquet",index=False)
M=M.merge(XA,left_on='a',right_index=True,how='left').merge(XL,left_on='a',right_index=True,how='left')
print('nonEUR done',time.time()-t0,flush=True)
M.to_parquet(W+"master_cs.parquet",index=False)
log['n_master']=len(M); log['elapsed_s']=time.time()-t0
json.dump(log,open(W+"s2_log.json","w"),indent=1,default=str); print(json.dumps(log,default=str))
