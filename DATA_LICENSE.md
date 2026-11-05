# Data licence

The derived resource (`resource/`), figure source data, results workbooks, analysis-plan documents and figures are released under the **Creative Commons Attribution 4.0 International licence (CC BY 4.0)**, https://creativecommons.org/licenses/by/4.0/.
Code is released under the MIT License (see `LICENSE`).

Please cite: Ni J, Zhang W, Zhang T. PAV-aware evidence tiers for UK Biobank Pharma Proteomics Project plasma cis-pQTLs: resource, code and analysis plans (v1.3.3). Zenodo; 2026. https://doi.org/10.5281/zenodo.23975177

## Third-party sources
Only derived per-credible-set statistics and annotations are released. **No third-party summary-statistic files are redistributed.** The only third-party table rows included are selected columns of Niu et al. 2025 Supplementary Table 4 (CC BY 4.0), listed below. Fields derived from the sources below remain subject to those sources' terms, and users must also cite the original studies.

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
| Generation Scotland MS pQTL summary statistics (Edinburgh DataShare doi:10.7488/ds/8011, CC BY 4.0) | `ms_gs_*`, `gs_*`, `ms_independent_combined` | CC BY 4.0 | Richmond et al. 2026 Nat Commun, doi:10.1038/s41467-026-76474-8 |
| HOLBAEK MS pQTL summary statistics, GWAS Catalog GCST90452968–GCST90453410 (all 170 studies used are CC0 1.0 per the GWAS Catalog API) | `ms_hb_*`, `hb_*`, `ms_independent_combined` | CC0 1.0 | Niu et al. 2025 Nat Genet, doi:10.1038/s41588-025-02089-2; Cerezo et al. 2025 |
| Niu et al. 2025 Supplementary Table 4 (pQTL list with peptide-evidence tiers; article CC BY 4.0) | `src/analysis/A4/niu_st4_tiers.tsv` (selected columns, redistributed with attribution; HOLBAEK sensitivity analysis) | CC BY 4.0 | Niu et al. 2025, doi:10.1038/s41588-025-02089-2 |
| Open Targets Platform 26.09 variant index | rsIDs in the lookup portal data (`docs/data/`) | Open Targets data, CC0 | Buniello et al. 2025 NAR, doi:10.1093/nar/gkae1128 |
