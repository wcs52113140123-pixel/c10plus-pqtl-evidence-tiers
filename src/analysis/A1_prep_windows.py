"""Addendum A1 prep: GRCh37 cis windows per UniProt accession and UKB-PPP cis CS variants lifted to GRCh37."""
import duckdb, pandas as pd, numpy as np
from pyliftover import LiftOver
B="../data_v1/db/parquet/"; W="./work/"
lo=LiftOver("../data_v1/raw/ref/hg38ToHg19.over.chain.gz")
con=duckdb.connect()
G=con.sql(f"select gene_id, symbol, chrom, gene_start, gene_end, uniprot_swissprot, uniprot_all from '{B}gene.parquet' where chrom in ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')").df()
def lift(ch,p):
    r=lo.convert_coordinate('chr'+ch,int(p)-1)
    if r and len(r)==1 and r[0][0]=='chr'+ch: return r[0][1]+1
    return None
rows=[]
for g in G.itertuples():
    s=lift(g.chrom,g.gene_start); e=lift(g.chrom,g.gene_end); L=g.gene_end-g.gene_start
    if s is None and e is None: continue
    if s is None: s=e-L
    if e is None: e=s+L
    a,b=min(s,e),max(s,e)
    w0=max(1,a-1_000_000-1000); w1=b+1_000_000+1000
    sw=list(g.uniprot_swissprot) if isinstance(g.uniprot_swissprot,(list,np.ndarray)) else []
    al=list(g.uniprot_all) if isinstance(g.uniprot_all,(list,np.ndarray)) else []
    for u in set(sw)|set(al): rows.append((u,g.gene_id,g.symbol,g.chrom,w0,w1, u in sw))
Wd=pd.DataFrame(rows,columns=['uniprot','gene_id','symbol','chrom','w0_37','w1_37','swissprot'])
# prefer swissprot mapping when an accession maps to several genes
Wd=Wd.sort_values(['uniprot','swissprot'],ascending=[True,False]).drop_duplicates('uniprot')
Wd.to_csv(W+"seer_windows_grch37.tsv",sep='\t',index=False)
U=con.sql(f"select gene_id from '{B}coverage_ukbppp_cis_gene_matrix.parquet'").df()
pd.Series(sorted(U.gene_id)).to_csv(W+"ukb_cis_genes.txt",index=False,header=False)
L=con.sql(f"""select l.studyLocusId a, c.gene_id ga, l.variant_id, l.pp, l.beta, l.is95, vc.chrom, vc.pos_b38, vc.ref, vc.alt
   from read_parquet('{B}qtl_cs_locus/*/*.parquet', hive_partitioning=true) l join '{B}qtl_cs.parquet' c using(studyLocusId) join '{B}variant_core.parquet' vc using(variant_id)
   where l.studyType='pqtl' and c.projectId='UKB_PPP_EUR' and not c.isTransQtl""").df()
L['pos37']=[lift(str(ch),p) for ch,p in zip(L.chrom,L.pos_b38)]
L.to_parquet(W+"ukb_cis_locus_grch37.parquet",index=False)
print(len(Wd), Wd.gene_id.nunique(), len(L), L.pos37.isna().sum())
