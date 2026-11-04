"""A4 prep: GS entry -> UKB cis gene mapping with GRCh37 windows; GRCh38 windows for all UKB cis genes (HOLBAEK)."""
import pandas as pd, numpy as np, duckdb, sys
W="./"; B="../data_v1/db/parquet/"
ent=pd.read_csv(W+'a4/gs_entries.tsv',sep='\t')
Wd=pd.read_csv(W+'work/seer_windows_grch37.tsv',sep='\t',dtype={'chrom':str})
res=pd.read_parquet(W+'a3/resource_v2.2.parquet',columns=['gene_id','symbol'])
ukb=set(res.gene_id)
u2g=Wd.set_index('uniprot')
rows=[]
for e in ent.itertuples():
    accs=[a.split('-')[0] for a in e.name.replace('.fastGWA.gz','').split('_')]
    gs=sorted({u2g.at[a,'gene_id'] for a in accs if a in u2g.index}&ukb)
    for g in gs:
        w=Wd[Wd.gene_id==g].iloc[0]
        rows.append(dict(name=e.name,lho=e.lho,csize=e.csize,usize=e.usize,crc=e.crc,gene_id=g,symbol=w.symbol,chrom=w.chrom,w0_37=int(w.w0_37),w1_37=int(w.w1_37),n_genes=len(gs),accs=';'.join(accs)))
G=pd.DataFrame(rows); G.to_csv(W+'a4/gs_entry_gene_windows.tsv',sep='\t',index=False)
print('entries mapped',G.name.nunique(),'genes',G.gene_id.nunique(),'multi',G[G.n_genes>1].name.nunique(), 'bytes',G.drop_duplicates('name').csize.sum()/1e9)
# GRCh38 windows for all UKB cis genes
gene=duckdb.sql(f"select gene_id,symbol,chrom,gene_start,gene_end,uniprot_swissprot,uniprot_all from '{B}gene.parquet'").df()
g38=gene[gene.gene_id.isin(ukb)].copy(); g38['w0_38']=(g38.gene_start-1_000_000-1000).clip(lower=1); g38['w1_38']=g38.gene_end+1_000_000+1000
g38['swiss']=g38.uniprot_swissprot.apply(lambda x:';'.join(x) if x is not None else ''); g38['allu']=g38.uniprot_all.apply(lambda x:';'.join(x) if x is not None else '')
g38[['gene_id','symbol','chrom','gene_start','gene_end','w0_38','w1_38','swiss','allu']].to_csv(W+'a4/ukb_cis_gene_windows_grch38.tsv',sep='\t',index=False)
print('grch38 windows',len(g38))
