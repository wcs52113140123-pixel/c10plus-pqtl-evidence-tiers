# Upstream data sources

| Source | Used for | Access / terms | Cite |
|---|---|---|---|
| UKB-PPP European-ancestry cis credible sets, Open Targets Platform 26.09 | discovery units, anchors (GTEx v10, eQTL Catalogue credible sets via Open Targets) | Open Targets data, CC0 | Sun et al. 2023 Nature, doi:10.1038/s41586-023-06592-6; Buniello et al. 2025 NAR, doi:10.1093/nar/gkae1128 |
| UKB-PPP non-European summary statistics, Synapse syn51364943 | cross-ancestry fields (`xa_*`) | registration required; the data page requires citation of Sun et al. 2023 | Sun et al. 2023, doi:10.1038/s41586-023-06592-6 |
| Sun et al. 2018 SomaScan v3 (INTERVAL) credible sets; eQTL Catalogue | `rep_sun2018`, `h4max_eqtl_eqtlcatalogue` | via Open Targets / eQTL Catalogue, CC BY 4.0 | Sun et al. 2018, doi:10.1038/s41586-018-0175-2; Kerimov et al. 2021, doi:10.1038/s41588-021-00924-w |
| GTEx v10 | `h4max_eqtl_gtex_v10`, wrong-gene control | public summary data | GTEx Consortium 2020, doi:10.1126/science.aaz1776 |
| deCODE SomaScan v4 (Eldjarn et al. 2023 supplementary tables) | `rep_decode_v4`, `dir_decode_v4` | publisher supplementary files | Eldjarn et al. 2023, doi:10.1038/s41586-023-06563-x |
| MS pQTL assessment (Suhre et al. 2025 supplementary tables) | `ms_primary`, `ms_mspa_*`, `ms_seer_discovery` | publisher supplementary files | Suhre et al. 2025, doi:10.1038/s41588-025-02413-w |
| Seer Proteograph summary statistics, GWAS Catalog GCST90570713–GCST90572692 | `ms_ext_*`, `seer_ext_*` | GWAS Catalog terms | Suhre et al. 2025; Cerezo et al. 2025, doi:10.1093/nar/gkae1070 |
| ARIC (Zhang et al. 2022 supplementary tables) | `aric_*`, `rep_aric_*`, `dir_aric_ea` | publisher supplementary files | Zhang et al. 2022, doi:10.1038/s41588-022-01051-w |
| Fenland (associations curated in the GWAS Catalog) | `fenland_has_cis`, `rep_fenland`, `dir_fenland` | GWAS Catalog terms | Pietzner et al. 2021 Science, doi:10.1126/science.abj1541; Cerezo et al. 2025 |
| Generation Scotland MS pQTLs (Edinburgh DataShare doi:10.7488/ds/8011; GWAS_sumstats.zip, 439 fastGWA entries, GRCh37) | `ms_gs_*`, `gs_*` | CC BY 4.0 | Richmond et al. 2026 Nat Commun, doi:10.1038/s41467-026-76474-8 |
| HOLBAEK MS pQTLs (GWAS Catalog GCST90452968–GCST90453410, harmonised GRCh38) | `ms_hb_*`, `hb_*` | CC0 1.0 (checked for all 170 studies used) | Niu et al. 2025 Nat Genet, doi:10.1038/s41588-025-02089-2 |
| Open Targets Platform 26.09 variant index (GraphQL `variant.rsIds`) | `rsIds` (portal data only; not a field of the resource TSV) | Open Targets data, CC0 | Buniello et al. 2025 NAR, doi:10.1093/nar/gkae1128 |

ARIC supplementary files were retrieved from the publisher (static-content.springer.com); Fenland associations from the GWAS Catalog REST API for PMID 34648354 (retrieved 2026-11-04). science.org and omicscience.org were not accessible programmatically. None of these files are redistributed here.

rsIDs shown and searchable on the lookup portal were taken from the Open Targets Platform 26.09 variant index (the same release as the credible sets); 14,076 of 14,078 lead variants have at least one rsID. They are added to the portal data file (`docs/data/tiers_v2.2.json`) only.

Generation Scotland windows were read from the DataShare ZIP with HTTP range requests (CRC-32 checked per entry); HOLBAEK windows were read with remote tabix (5 studies verified against full-file downloads and GWAS Catalog md5). Only per-credible-set derived fields are released; the extracted windows are not redistributed.
