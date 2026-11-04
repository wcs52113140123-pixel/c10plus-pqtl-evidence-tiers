
import duckdb, pandas as pd
B="../data_v1/db/parquet/"
con=duckdb.connect()
print(con.sql(f"select layer, cis_trans_src, is_cis_1mb, count(*) from '{B}ext_pqtl_layers.parquet' where layer like 'decode_somascan%' or layer like 'seer%' group by all order by 1").df().to_string())
print(con.sql(f"select cs_members_src from '{B}ext_pqtl_layers.parquet' where layer='decode_somascan_ST33_cs' and cs_members_src is not null limit 2").df().to_string(max_colwidth=300))
print(con.sql(f"select unnest(alleleFrequencies).populationName p, count(*) from (select * from '{B}variant_core.parquet' limit 20000) group by 1").df().to_string())
print(con.sql(f"select studyType, projectId, count(*), sum((l.beta is null)::int) nb, sum((l.logBF is null)::int) nl from read_parquet('{B}qtl_cs_locus/*/*.parquet', hive_partitioning=true) l join '{B}qtl_cs.parquet' c using(studyLocusId) where projectId in ('UKB_PPP_EUR','Sun_2018','GTEx-v10') group by all").df().to_string())
print(con.sql(f"select projectId, count(*) , count(distinct gene_id), sum((assay_id is null)::int) from '{B}pqtl_study.parquet' group by 1").df())
print(con.sql(f"select studyId, gene_id, assay_id, platform from '{B}pqtl_study.parquet' where projectId='UKB_PPP_EUR' limit 3").df())
print(con.sql(f"select gene_id, count(distinct oid) n from read_parquet('{B}ukbppp_nonEUR_cis_sumstats/*/*.parquet', hive_partitioning=true) group by 1 having n>1 limit 5").df())
