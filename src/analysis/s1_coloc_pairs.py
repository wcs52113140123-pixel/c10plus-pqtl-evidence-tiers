"""Stage 1: QTL x QTL colocalisation sufficient statistics for every UKB-PPP EUR cis credible set.
For each pair (a = UKB cis CS, b = comparator CS) sharing >=1 variant, compute the log-sum-exp aggregates
needed for coloc single-effect PP.H4 under (i) the pilot 'union' convention (variant absent from one CS -> lbf 0)
and (ii) the intersection (coloc.bf_bf) convention. Pairs sharing no variant have PP.H4 ~ 0 and are not stored.
Comparators: Sun_2018 cis pQTL CS (cognate gene), all eQTL/sc-eQTL/sQTL/tuQTL CS of the cognate gene,
GTEx-v10 eQTL CS of non-cognate genes (wrong-gene control), UKB-PPP trans CS (hotspot sensitivity)."""
import duckdb, time, numpy as np, pandas as pd, os
t0=time.time()
B="../data_v1/db/parquet/"
OUT="./work/"; os.makedirs(OUT,exist_ok=True)
con=duckdb.connect(); con.sql("SET threads=96; SET memory_limit='180GB'")
con.sql(f"create view cs as select * from '{B}qtl_cs.parquet'")
con.sql(f"create view loc as select * from read_parquet('{B}qtl_cs_locus/*/*.parquet', hive_partitioning=true)")
con.sql("create table ukb as select studyLocusId a, gene_id ga from cs where projectId='UKB_PPP_EUR' and not isTransQtl")
con.sql("""create table cmeta as select studyLocusId b, studyType, projectId, gene_id gb, biosampleId, studyId,
           (projectId='UKB_PPP_EUR' and isTransQtl) as ukbtrans from cs
           where studyType<>'pqtl' or (projectId='Sun_2018' and not isTransQtl) or (projectId='UKB_PPP_EUR' and isTransQtl)""")
con.sql("create table la as select l.studyLocusId a, l.variant_id v, l.logBF l1, l.pp pp1, l.beta beta1, l.is95 is95a from loc l join ukb u on l.studyLocusId=u.a where l.studyType='pqtl'")
print('ukb cs',con.sql('select count(*) from ukb').fetchone(),'la rows',con.sql('select count(*), sum((l1 is null)::int) from la').fetchone(),flush=True)
con.sql("create table vset as select distinct v from la")
con.sql("""create table lb as select l.studyLocusId b, l.variant_id v, l.logBF l2, l.pp pp2, l.beta beta2 from loc l
           join vset s on l.variant_id=s.v join cmeta m on l.studyLocusId=m.b""")
print('lb rows (shared-variant rows)',con.sql('select count(*) from lb').fetchone(),time.time()-t0,flush=True)
# log-sum-exp in two passes
con.sql("""create table j as select la.a, lb.b, la.v, l1, l2, pp1, pp2, beta1, beta2 from la join lb using(v)""")
con.sql("""create table pairs as with mx as (select a,b,max(l1+l2) m12,max(l1) m1,max(l2) m2,count(*) k from j group by a,b)
 select j.a,j.b,mx.k, mx.m12+ln(sum(exp(l1+l2-mx.m12))) Lsh12, mx.m1+ln(sum(exp(l1-mx.m1))) Lsh1, mx.m2+ln(sum(exp(l2-mx.m2))) Lsh2,
   arg_max(j.v, pp1*pp2) vtop, arg_max(beta1, pp1*pp2) b1top, arg_max(beta2, pp1*pp2) b2top
 from j join mx using(a,b) group by j.a,j.b,mx.k,mx.m12,mx.m1,mx.m2""")
con.sql("""create table pairs2 as select p.*, m.studyType, m.projectId, m.gb, m.biosampleId, m.studyId, m.ukbtrans, u.ga from pairs p join cmeta m using(b) join ukb u using(a)
           where m.gb=u.ga or (m.projectId='GTEx-v10' and m.studyType='eqtl') or m.ukbtrans""")
print('pairs',con.sql('select count(*) from pairs2').fetchone(),time.time()-t0,flush=True)
# per-CS totals over all CS variants
con.sql("create table bneed as select distinct b from pairs2")
con.sql("""create table tb as with x as (select l.studyLocusId b, l.logBF l2 from loc l join bneed n on l.studyLocusId=n.b)
   , mx as (select b, max(l2) mb, count(*) nb from x group by b)
   select x.b, mx.nb, mx.mb+ln(sum(exp(l2-mx.mb))) Sb from x join mx using(b) group by x.b, mx.nb, mx.mb""")
con.sql("""create table ta as with mx as (select a, max(l1) ma, count(*) na from la group by a)
   select la.a, mx.na, mx.ma+ln(sum(exp(l1-mx.ma))) Sa from la join mx using(a) group by la.a, mx.na, mx.ma""")
df=con.sql("select p.*, ta.na, ta.Sa, tb.nb, tb.Sb from pairs2 p join ta using(a) join tb using(b)").df()
df.to_parquet(OUT+"pairs_stats.parquet",index=False)
con.sql("select a, ga from ukb").df().to_parquet(OUT+"ukb_cis_cs.parquet",index=False)
print('saved',len(df),df.studyType.value_counts().to_dict(),time.time()-t0,flush=True)
