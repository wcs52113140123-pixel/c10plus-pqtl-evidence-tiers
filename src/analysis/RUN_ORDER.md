# C10+ full-scale analysis — run order (a local Linux workstation, venv c10_pqtl/venv; seed 20261103)
Data: ../data_v1/db/parquet ; work dir analysis_v2/
1. s1_coloc_pairs.py   QTLxQTL coloc sufficient statistics for every UKB-PPP cis CS (pairs sharing >=1 variant)
2. s2_build_master.py  per-CS master table (anchors, PAV, Sun/deCODE replication, MSPA/Seer MS truth, non-EUR meta); needs work/suhre2025_MSPA.xlsx
3. h1_h2.py            H1 (Sun 2018) and H2a/H2b/H2c (deCODE v4)
4. h3.py               H3 mass-spec truth
5. h4.py               H4 grade models, tiers, final scores
6. h5.py               H5 cross-ancestry
7. sensitivity.py      pre-specified sensitivity analyses
8. build_resource.py   resource table + cardiometabolic application
Post hoc (beta-consistency flag, slope robustness) and figures were run in the interactive kernel; their source data are the *_source.csv files.
Software: python 3.10, duckdb 1.5.6, pandas 2.3.3, numpy 2.2.6, scipy 1.15.3, statsmodels 0.14.6, scikit-learn 1.7.2.

## Addendum A1 (Seer full)
9. A1_prep_windows.py (a local Linux workstation): GRCh37 windows per UniProt + UKB cis CS variants lifted to GRCh37 (hg38ToHg19 chain)
10. A1_seer_worker.sh / A1_seer_worker2.sh (a rented cloud server, sandbox): streaming cis-window extraction (v2 = early stop on sorted files)
11. A1_validate.py: structural validation of extracted windows
12. A1_ms_ext.py offset -> offset_choice.json; A1_ms_ext.py full -> labels, H3i-ext/H3ii-ext/H4-ext

## Addendum A2 (ARIC + Fenland)
13. A2_prep_A2.py: ARIC signal tables (ST3/ST6/ST10), Fenland name->gene mapping (SomaScan names, then UniProt Swiss-Prot names; uniprot_human_sp_names.tsv.gz from rest.uniprot.org)
14. A2_a2_analysis.py: H2a/H2b per cohort, effect direction, DL meta, per-CS fields

## Paths
Scripts use paths relative to the analysis working directory (`./work`, `./results`, `./seer_ext`) and expect the raw database under `../data_v1/` (DuckDB/parquet tables built from the sources listed in docs/DATA_SOURCES.md). The figure script reads its inputs from `inputs/` (results workbook, resource table and ROC source CSV, as listed at the top of the script).
