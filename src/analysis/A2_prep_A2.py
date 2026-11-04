"""A2 prep: ARIC signal tables and Fenland trait->gene mapping (counts and mapping rates only; no outcome computed)."""
import pandas as pd, numpy as np, re, json, gzip, duckdb
MAN="../data_v1/raw/manual/"; B="../data_v1/db/parquet/"
W="./work/"
con=duckdb.connect()
G=con.sql(f"select gene_id, symbol, chrom, gene_start, gene_end, uniprot_swissprot from '{B}gene.parquet'").df()
sym2g=G.groupby('symbol').gene_id.apply(list).to_dict()
u2g={}
for g,us in zip(G.gene_id,G.uniprot_swissprot):
    if isinstance(us,(list,np.ndarray)):
        for u in us: u2g.setdefault(u,set()).add(g)
X=pd.ExcelFile(MAN+"aric/MOESM4_ESM.xlsx")
def sheet(name,hdr):
    d=pd.read_excel(X,sheet_name=name,header=hdr); return d
log={}
def gene_of(uni,sym):
    gs=set()
    for u in re.split(r'[|,; ]+',str(uni)):
        gs|=u2g.get(u,set())
    if len(gs)==1: return list(gs)[0]
    s=sym2g.get(str(sym),[])
    return s[0] if len(s)==1 else None
out={}
for anc,s3,s6,s10 in [('EA','ST3.1--significant SOMAmer (EA)','ST6.1 -- cond-indep pQTL (EA)','ST10.1-- fine mapping (EA)'),('AA','ST3.2--significant SOMAmer (AA)','ST6.2 -- cond-indep pQTL (AA)','ST10.2-- fine mapping (AA)')]:
    a=sheet(s3,1).iloc[1:]; a.columns=[str(c) for c in a.columns]
    a=a.rename(columns={a.columns[0]:'somamer',a.columns[1]:'target',a.columns[2]:'fullname',a.columns[3]:'uniprot',a.columns[4]:'symbol',a.columns[5]:'chr',a.columns[8]:'rsid',a.columns[9]:'pos',a.columns[10]:'allele',a.columns[11]:'af',a.columns[12]:'beta',a.columns[13]:'p'})
    a=a[a.somamer.astype(str).str.startswith('SeqId')]
    b=sheet(s6,1); b.columns=['somamer','target','fullname','uniprot','symbol','chr','tss','sentinel','sentinel_pos','rsid','pos','allele','af','rank','ld','beta_joint','p_joint','beta','p']
    b=b[b.somamer.astype(str).str.startswith('SeqId')]
    f=sheet(s10,1); f=f.iloc[:,:9]; f.columns=['somamer','symbol','chr','rsid','pos','pip','p','ld_sentinel','dist_kb']
    f=f[f.somamer.astype(str).str.startswith('SeqId')]
    som=pd.concat([a[['somamer','uniprot','symbol','target','fullname']],b[['somamer','uniprot','symbol','target','fullname']]]).drop_duplicates('somamer')
    som['gene_id']=[gene_of(u,s) for u,s in zip(som.uniprot,som.symbol)]
    smap=som.set_index('somamer').gene_id.to_dict()
    S=pd.concat([a.assign(src='sentinel',pip=np.nan)[['somamer','chr','rsid','pos','allele','af','beta','p','src','pip']],
                 b.assign(src='conditional',pip=np.nan)[['somamer','chr','rsid','pos','allele','af','beta','p','src','pip']],
                 f.assign(src='finemap',allele=np.nan,af=np.nan,beta=np.nan)[['somamer','chr','rsid','pos','allele','af','beta','p','src','pip']]])
    S['gene_id']=S.somamer.map(smap); S['chr']=S.chr.astype(str); S['pos']=pd.to_numeric(S.pos,errors='coerce')
    S.to_csv(W+f"a2_aric_signals_{anc}.tsv",sep='\t',index=False); som.to_csv(W+f"a2_aric_somamers_{anc}.tsv",sep='\t',index=False)
    log[f'aric_{anc}']=dict(somamers_with_cis=int(a.somamer.nunique()),somamers_mapped=int(som.gene_id.notna().sum()),somamers_total_listed=int(len(som)),genes=int(som.gene_id.nunique()),
        n_sentinel=int((S.src=='sentinel').sum()),n_conditional=int((S.src=='conditional').sum()),n_finemap=int((S.src=='finemap').sum()))
    out[anc]=som
# ---- Fenland mapping ----
def norm(s): 
    s=str(s).replace('\u2019',"'").replace('\u2018',"'").strip(); s=re.sub(r'\s+levels?$','',s,flags=re.I); return re.sub(r'\s+',' ',s).strip().lower()
def strip_par(s): return re.sub(r'\s*\([^()]*\)\s*$','',s).strip()
D1={}
for anc,som in out.items():
    for nm,g in list(zip(som.fullname,som.gene_id))+list(zip(som.target,som.gene_id)):
        if isinstance(nm,str) and g: D1.setdefault(norm(nm),set()).add(g)
U=pd.read_csv(MAN+"uniprot_human_sp_names.tsv.gz",sep='\t')
D2={}
for acc,pn in zip(U.Entry,U['Protein names']):
    gs=u2g.get(acc,set())
    if len(gs)!=1 or not isinstance(pn,str): continue
    g=list(gs)[0]
    pn=re.sub(r'\[.*?\]','',pn)
    rec=pn.split(' (')[0]; alts=re.findall(r'\(([^()]*)\)',pn)
    for nm in [rec]+[a for a in alts if not a.startswith('EC ')]:
        D2.setdefault(norm(nm),set()).add(g)
St=pd.read_csv(MAN+"fenland/fenland_gwascat_studies.tsv",sep='\t')
def mapname(t):
    t0=norm(t); t1=strip_par(t0)
    complex_=bool(re.search(r':|heterodimer|heterotrimer|heterotetramer|complex',t0))
    for D,lab in [(D1,'somascan_name'),(D2,'uniprot_name')]:
        for k,kl in [(t0,'exact'),(t1,'strip_paren')]:
            if k in D:
                gs=D[k]
                if len(gs)==1: return list(gs)[0], f'{lab}_{kl}', complex_
                return None, f'ambiguous_{lab}', complex_
    return None, 'unmapped', complex_
M=[mapname(t) for t in St.disease_trait]
St['gene_id']=[m[0] for m in M]; St['map_method']=[m[1] for m in M]; St['complex']=[m[2] for m in M]
St.loc[St.complex & ~St.map_method.str.startswith('somascan'),'gene_id']=None
St.to_csv(W+"a2_fenland_study_gene_map.tsv",sep='\t',index=False)
A=pd.read_csv(MAN+"fenland/fenland_gwascat_associations.tsv",sep='\t')
A=A.merge(St[['accession_id','gene_id','map_method']],on='accession_id',how='left')
loc=A.location_grch38.astype(str).str.split(':',expand=True); A['chr']=loc[0]; A['pos']=pd.to_numeric(loc[1],errors='coerce')
A=A.merge(G[['gene_id','chrom','gene_start','gene_end']],on='gene_id',how='left')
A['cis']=(A.chr==A.chrom)&(A.pos>=A.gene_start-1e6)&(A.pos<=A.gene_end+1e6)
A.to_csv(W+"a2_fenland_assoc_mapped.tsv",sep='\t',index=False)
ukb=set(open(W+'ukb_cis_genes.txt').read().split())
log['fenland']=dict(studies=int(len(St)),mapped=int(St.gene_id.notna().sum()),map_rate=float(St.gene_id.notna().mean()),by_method=St.map_method.value_counts().to_dict(),complex_excluded=int((St.complex & ~St.map_method.str.startswith('somascan')).sum()),
   genes=int(St.gene_id.nunique()),assoc=int(len(A)),assoc_mapped=int(A.gene_id.notna().sum()),cis_assoc=int(A.cis.sum()),genes_with_cis=int(A[A.cis].gene_id.nunique()),ukb_cis_genes_with_fenland_cis=int(len(set(A[A.cis].gene_id)&ukb)))
for anc in ['EA','AA']:
    S=pd.read_csv(W+f"a2_aric_signals_{anc}.tsv",sep='\t'); log[f'aric_{anc}']['ukb_cis_genes_with_aric_cis']=int(len(set(S.gene_id.dropna())&ukb))
json.dump(log,open(W+"a2_prep_log.json","w"),indent=1,default=str); print(json.dumps(log,indent=1,default=str))
