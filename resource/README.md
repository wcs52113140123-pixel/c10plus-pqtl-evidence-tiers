# C10+ PAV-aware evidence tiers for UKB-PPP cis-pQTL credible sets (v2.3)

**Files**: `C10plus_cis_pqtl_robustness_grades_v2.3.tsv.gz` (SHA-256 `067b7f7c44767b1befa429e5721dafac8e7d5ecd8a88473ffd2e7ce550b3c8be`) and `.parquet`
**Rows**: 14,472 = all Open Targets 26.09 UKB-PPP European-ancestry cis credible sets (SuSiE-inf) from the UKB-PPP Olink Explore 3072 study, 1,879 genes; 101 columns. Replication layers: SomaScan v3 (Sun et al. 2018, INTERVAL), SomaScan v4 (Eldjarn et al. 2023, deCODE; assay availability inferred from the published signal tables) and, in addendum A2, SomaScan v4 ARIC (Zhang et al. 2022; European and African ancestry) and Fenland (Pietzner et al. 2021; GWAS Catalog lead associations).
**v2.3 vs v2.2**: the 76 v2.2 columns are unchanged value for value; 25 addendum A4 columns were appended with mass-spectrometry labels from two cohorts that are independent of the Seer study: Generation Scotland (Richmond et al. 2026; serum Scanning SWATH; N = 15,035; primary) and HOLBAEK (Niu et al. 2025; plasma Orbitrap Astral DIA; N = 1,909; secondary). Tiers and concordance are unchanged. **v2.2 vs v2.1**: the 67 v2/v2.1 columns are unchanged value for value; nine addendum A2 columns were appended (see the end of the field table). Tiers, concordance and MS labels do not use the A2 layers. (v2.1 vs v2: README and example scripts only.)
**Definitions**: frozen pre-registration v2.0 (SHA-256 126f515f4734baef319629584521bc492e4a70730067f293c76964c2b62a89b0) addendum A1 (SHA-256 6d54bd5478c43184e349523b95c0406ded7bab5e106823e0d7537ad85b2035e2) addendum A2 (SHA-256 08df92f09ddaa5e82e54e87e134cab325dc6d02224b00456dcde6599522cbdb3), addendum A3 (SHA-256 4a5170ac229e012b253a36046f7163144e0669ec1fec8f66b4e77d81279c3971) and addendum A4 (SHA-256 cc2da0441cba82e1248ff6c328ac692cc5043fa0b6c133a5eb30f74c01935fd7).
**Licence**: CC BY 4.0 · **Citation**: Ni J, Zhang W, Zhang T. Protein-altering variants inflate cross-platform replication of plasma protein QTLs (manuscript in submission). Data: https://doi.org/10.5281/zenodo.23975177 (version 1.3.3; all versions https://doi.org/10.5281/zenodo.23938290); code: GitHub tag v1.3.3. Portal: https://wcs52113140123-pixel.github.io/c10plus-pqtl-evidence-tiers/ (rsIDs shown on the portal come from Open Targets 26.09).

## How to use
- `tier` is the pre-registered rule-based grade and does not depend on model training:
  - **T1** eQTL-anchored (coloc PP.H4 ≥ 0.8 with cognate GTEx v10 or eQTL Catalogue eQTL), no PAV in the 95% set — 2,539 CS
  - **T2** anchored with PAV, or isoform-anchored (sQTL/tuQTL) without PAV — 682 CS
  - **T3** no anchor, no PAV — 10,046 CS
  - **T4** no eQTL anchor, PAV in the 95% set, regardless of platform concordance — 1,205 CS
- **T4 does not mean false positive.** It means affinity measurements alone cannot separate an abundance effect from an epitope effect (functional coding variants in PCSK9, APOE and LPL are T4). Seek MS, functional or orthogonal biomarker evidence.
- `flag_concordant_unanchored_PAV` (339 CS, 254 genes): T4 signals that replicated on SomaScan; concordance is least informative here.
- `platform_concordance` is informational: PAVs inflated binary Olink–SomaScan replication (power-matched OR 2.74), so concordance is not orthogonal validation for PAV-containing signals. In general, concordance was associated with MS support in the two independent MS cohorts (combined OR 2.92) but only weakly in the Seer layers.
- `p_ms_support_model` has poor external validity (held-out AUC 0.68; below-chance transfer to other labels); do not use it alone for filtering.
- `flag_beta_inconsistent` (204 CS, post hoc): OT lead β > 5× the value expected from P and MAF; exclude before effect-size analyses.
- Colocalisation uses per-variant logBF within OT 99% sets (approximation to full-region coloc). MAF is gnomAD NFE (UKB-PPP eaf unavailable).
- MS columns: `ms_primary`, `ms_mspa_peptide`, `ms_seer_discovery` (Suhre et al. 2025 tables) and the v2 fields `ms_ext_primary`, `ms_ext_nominal`, `ms_ext_cswide`, `seer_ext_*` (addendum A1: Seer full summary statistics, 564 studies; 167 supported / 200 refuted / 1,344 indeterminate). Both Seer label sets come from one MS study (agreement κ = 0.89) and are not independent. v2.3 adds `ms_gs_*`/`gs_*` (Generation Scotland) and `ms_hb_*`/`hb_*` (HOLBAEK) labels computed with the A1 rule; these two layers agree with each other (κ = 0.57) but hardly with the Seer labels (κ −0.02 to 0.11). eQTL anchoring did not predict MS support in any of the four layers; PAV-containing CS were less often MS-supported in all four (independent cohorts combined OR 0.61, P = 0.056; secondary analysis). No MS layer includes chromosome X (500 CS carry no MS label).

## Interpreting tiers (read before filtering)
Tier numbers code **evidence patterns, not an ordinal validity scale**. Observed rates by tier (T1/T2/T3/T4):

| Outcome | T1 | T2 | T3 | T4 | Trend P |
|---|---|---|---|---|---|
| MS supported, MSPA | 49/92 (53%) | 14/35 (40%) | 16/24 (67%) | 4/18 (22%) | 0.17 |
| MS supported, Seer full (A1) | 86/179 (48%) | 24/61 (39%) | 49/103 (48%) | 8/24 (33%) | 0.39 |
| MS supported, Generation Scotland (A4) | 50/78 (64%) | 18/31 (58%) | 221/303 (73%) | 24/38 (63%) | 0.25 |
| MS supported, HOLBAEK (A4) | 123/147 (84%) | 37/47 (79%) | 164/185 (89%) | 19/24 (79%) | 0.42 |
| Cross-ancestry replicated | 1,681/2,387 (70%) | 433/616 (70%) | 2,078/4,594 (45%) | 345/468 (74%) | non-monotonic |

- MS denominators are small and selected: only 0.2% (T3) to 5.1% (T2) of credible sets are MSPA-evaluable.
- T1 means "eQTL colocalisation and no PAV". This is a mechanistic argument against a binding-only explanation, not evidence from MS that the signal is true.
- MS `refuted` labels (both layers) assume equal standardized effects on Olink and Seer. In the A1 layer, observed Seer |z| was a median 0.50 of expectation among well-powered credible sets, and refutations nearly vanish if expected effects are reduced by 10–20%. Read `refuted` as "not detected despite nominal power".
- Peptide-level MS artefacts are possible for PAV signals (for example, the APOE ε2/ε3/ε4 residues alter tryptic cleavage). `ms_pav_peptide_flag` covers MSPA only; no equivalent exists for A1.
- Ties: several genes (for example, F11, APOE, KLKB1, LGALS3, LPA) have two or more credible sets at the capped −log10 P (≈322). The examples break ties by larger |lead_beta| and print a warning. Inspect tied credible sets individually.

## Versioning and reproducibility
- Anchors depend on GTEx v10 and the eQTL Catalogue as distributed in Open Targets 26.09, and replication depends on the published SomaScan tables. Later releases may change anchors and replication calls. Any re-grading will be released as a new version with its own checksums. Tier rules will not change without a new pre-registration.
- Analysis code (with run order) and figure code are provided with the Zenodo archive.
- A4 caveat: Generation Scotland and HOLBAEK labels use the equal-standardized-effect power rule of A1; read `refuted` as "not detected despite nominal power". MS quantification is itself sequence-dependent.
- A2 caveat: deCODE, ARIC and Fenland all use SomaScan v4, so agreement among them is cross-cohort reproducibility, not platform orthogonality. Do not count `rep_aric_ea`/`rep_fenland` as independent validation of an epitope-sensitive signal.

## Examples
(`FIELD_DICTIONARY_zh.md` is a Chinese field dictionary and `FIELD_DICTIONARY_en.md` its English translation; for field definitions, this README is authoritative.)
- Python: `python example_select_instruments.py PCSK9`
- R: `Rscript example_select_instruments.R PCSK9`
Both list all cis CS for a gene (ties at the capped P broken by |β|, with a warning), report the strongest Tier 1 credible set as a candidate sensitivity instrument and summarise flags (tested with pandas 2.3 / R 4.x readr+dplyr).

## Fields
| field | description |
|---|---|
| `studyLocusId` | Open Targets 26.09 credible-set ID (UKB_PPP_EUR, isTransQtl=false) |
| `studyId` | OT study ID (UKB-PPP Olink assay) |
| `gene_id` | Ensembl gene (cognate) |
| `symbol` | HGNC symbol |
| `olink_oid` | Olink assay OID |
| `chrom` | chromosome |
| `lead_variant_id` | lead variant chrom_pos_ref_alt (GRCh38) |
| `pos_b38` | lead position |
| `lead_beta` | OT lead beta (alt allele; NPX scale) |
| `lead_mlog10p` | −log10 P of lead |
| `lead_maf_gnomad_nfe` | lead MAF from gnomAD NFE (OT variant index; UKB eaf unavailable) |
| `n_variants_cs95` | variants in 95% credible set |
| `purity_min_r2` | SuSiE purity (min r²) |
| `dist_lead_to_tss` | |lead − cognate TSS| (bp) |
| `h4max_eqtl_gtex_v10` | max coloc PP.H4 vs cognate GTEx v10 eQTL CS (any tissue; SuSiE logBF; p1=p2=1e-4, p12=1e-5) |
| `h4max_eqtl_eqtlcatalogue` | max PP.H4 vs cognate eQTL Catalogue eQTL/sc-eQTL CS |
| `h4max_sqtl_tuqtl` | max PP.H4 vs cognate sQTL/tuQTL CS (GTEx + eQTL Catalogue) |
| `h4max_wronggene_eqtl_gtex` | max PP.H4 vs non-cognate GTEx eQTL CS (negative control) |
| `n_anchor_eqtl_datasets` | number of eQTL studies with PP.H4≥0.8 |
| `anchor_eqtl` | 1 if GTEx or eQTL Catalogue eQTL PP.H4≥0.8 |
| `anchor_eqtl_gtex` | GTEx eQTL anchor |
| `anchor_eqtl_nongtex` | eQTL Catalogue anchor |
| `anchor_isoform` | sQTL/tuQTL PP.H4≥0.8 |
| `pav_cognate_cs95_vep116_alltx` | 95% CS contains cognate-gene protein-altering variant (VEP 116, any transcript) — primary PAV |
| `pav_pip_sum` | summed PIP of cognate PAVs |
| `pav_cognate_cs95_canonical` | PAV on canonical transcript |
| `pav_lead_variant` | lead variant is cognate PAV |
| `somascan_v3_sun2018_assayed` | gene measured by SomaScan v3 (Sun 2018) |
| `h4max_sun2018` | max PP.H4 vs Sun 2018 cis CS |
| `rep_sun2018` | Sun 2018 replication (PP.H4≥0.8); NA if not assayed |
| `dir_sun2018` | effect direction at top shared variant (same/opposite/undetermined/NA) |
| `somascan_v4_decode_assayed` | gene present in Eldjarn 2023 deCODE SomaScan v4 tables |
| `rep_decode_v4` | deCODE cis signal (lead or LD-class member) overlaps UKB 95% CS; NA if not assayed |
| `dir_decode_v4` | direction at deCODE lead variant |
| `platform_concordance` | concordant (replicated on ≥1 SomaScan layer) / discordant (assayed, not replicated) / na (not assayed) |
| `ms_mspa_source` | Suhre 2025 MSPA table matched (ST10 UKB-PPP / ST9 deCODE) |
| `ms_primary` | supported (Bonferroni MS replication) / refuted (power>99%, not replicated) / indeterminate_underpowered / not_assessed |
| `ms_mspa_peptide` | MSPA>0.8 supported / <0.2 refuted / intermediate / not_assessed |
| `ms_pav_peptide_flag` | MSPA row lists a PAV-containing peptide |
| `ms_seer_discovery` | Seer MS cis-pQTL (ST1–ST3) allele in UKB 99% CS: overlap/no_overlap/not_assessed |
| `xa_n_ancestries` | non-EUR UKB-PPP ancestries contributing |
| `xa_meta_beta_aligned` | fixed-effect meta beta aligned to EUR direction |
| `xa_meta_se` | meta SE |
| `xa_meta_p` | meta P |
| `xa_replicated` | aligned beta>0 and P<0.05 |
| `xa_same_direction` | aligned beta>0 |
| `in_hla` | lead in chr6:28.51–33.48 Mb (none in OT set) |
| `n_other_genes_trans_cs_shared` | other proteins whose UKB-PPP trans CS share a variant (hotspot ≥5) |
| `tier` | pre-registered rule tier 1–4 |
| `tier_label` | T1 anchored no PAV; T2 anchored+PAV or isoform-only; T3 no anchor no PAV; T4 no anchor + PAV (regardless of platform concordance) |
| `p_ms_support_model` | L2 logistic (all chromosomes, MS target) predicted probability — low external validity, see H4 |
| `p_xa_rep_model` | L2 logistic (xa target) predicted probability |
| `model_confidence` | tertile of p_ms_support_model |
| `flag_concordant_unanchored_PAV` | platform-concordant but unanchored and PAV-containing (key caution class) |
| `lead_beta_to_expected_ratio` | |beta| / (z / sqrt(2pqN), N=34,000) |
| `flag_beta_inconsistent` | ratio>5 (OT beta inconsistent with P/MAF; 204 CS) — post hoc |
| `ms_ext_primary` | A1 primary label: supported (index variant Seer P<0.05/1,779, concordant direction) / refuted (power ≥0.99, not supported) / indeterminate / not_assessed (v2) |
| `ms_ext_nominal` | as above at P<0.05 with power at α=0.05 (v2) |
| `ms_ext_cswide` | supported if any evaluable 95%-CS variant has P<α/n_eval with concordant direction (v2) |
| `seer_ext_opposite_sig` | a significant opposite-direction (CS, study) pair exists (v2) |
| `seer_ext_studies` | GCST studies evaluated (v2) |
| `seer_ext_index_variant` | index variant (highest-PIP evaluable 95%-CS variant) in the best-powered study (v2) |
| `seer_ext_beta_aligned` | Seer effect aligned to the OT alt allele (inverse-normal scale) (v2) |
| `seer_ext_se` | Seer SE (v2) |
| `seer_ext_p` | Seer P (v2) |
| `seer_ext_power` | power at the Bonferroni threshold (v2) |
| `seer_ext_zexp` | expected z = |z_UKB|·√(N_Seer/N_UKB) (v2) |
| `aric_ea_has_cis` | 1 if the cognate gene has ≥1 published ARIC European-ancestry cis signal (SomaScan v4; A2 analysis set) (v2.2) |
| `rep_aric_ea` | 1 if any ARIC EA sentinel, conditionally independent or fine-mapped cis variant of the cognate gene lies in the 95% CS (exact position; no LD proxies); NA if `aric_ea_has_cis` = 0 (v2.2) |
| `rep_aric_ea_strict` | as `rep_aric_ea` but only ARIC fine-mapped variants with PIP ≥ 0.5 (v2.2) |
| `dir_aric_ea` | effect direction of the matched ARIC EA variant vs UKB-PPP after allele harmonisation: same / opposite / undetermined / NA (v2.2) |
| `aric_aa_has_cis` | as `aric_ea_has_cis` for ARIC African ancestry (secondary) (v2.2) |
| `rep_aric_aa` | as `rep_aric_ea` for ARIC African ancestry (secondary) (v2.2) |
| `fenland_has_cis` | 1 if the cognate gene has ≥1 GWAS Catalog Fenland cis lead association after name-based gene mapping (93.0% of studies mapped) (v2.2) |
| `rep_fenland` | 1 if any Fenland cis lead variant of the cognate gene lies in the 95% CS (exact position; no LD proxies); NA if `fenland_has_cis` = 0 (v2.2) |
| `dir_fenland` | effect direction of the matched Fenland lead vs UKB-PPP: same / opposite / NA (v2.2) |

| `ms_gs_primary` | Generation Scotland primary label (A1 rule; K = 1,103 pairs, α = 4.5×10⁻⁵): supported / refuted / indeterminate / not_assessed (v2.3) |
| `ms_gs_nominal` | as above at P < 0.05 with power at α = 0.05 (v2.3) |
| `ms_gs_cswide` | supported if any evaluable 95%-CS variant has P < α/n_eval with concordant direction (v2.3) |
| `gs_opposite_sig` | a significant opposite-direction (CS, entry) pair exists (v2.3) |
| `gs_entries` | Generation Scotland DataShare entries (UniProt-group fastGWA files) evaluated (v2.3) |
| `gs_index_variant` | index variant (highest-PIP evaluable 95%-CS variant) (v2.3) |
| `gs_beta_aligned` | Generation Scotland effect aligned to the OT alt allele (v2.3) |
| `gs_se` | Generation Scotland SE (v2.3) |
| `gs_p` | Generation Scotland P (v2.3) |
| `gs_power` | power at the Bonferroni threshold (v2.3) |
| `gs_zexp` | expected z = |z_UKB|·√(N_GS/N_UKB) (v2.3) |
| `ms_hb_primary` | HOLBAEK primary label (A1 rule; K = 720 pairs, α = 6.9×10⁻⁵): supported / refuted / indeterminate / not_assessed / not_evaluable (no public HOLBAEK summary statistics for the protein) (v2.3) |
| `ms_hb_nominal` | as above at P < 0.05 (v2.3) |
| `ms_hb_cswide` | CS-wide HOLBAEK label (v2.3) |
| `hb_opposite_sig` | a significant opposite-direction (CS, study) pair exists (v2.3) |
| `hb_entries` | HOLBAEK GWAS Catalog studies (GCST90452968–GCST90453410) evaluated (v2.3) |
| `hb_index_variant` | index variant (v2.3) |
| `hb_beta_aligned` | HOLBAEK effect aligned to the OT alt allele (v2.3) |
| `hb_se` | HOLBAEK SE (v2.3) |
| `hb_p` | HOLBAEK P (v2.3) |
| `hb_power` | power at the Bonferroni threshold (v2.3) |
| `hb_zexp` | expected z = |z_UKB|·√(N_HB/N_UKB) (v2.3) |
| `hb_evaluable` | True if the cognate protein has a HOLBAEK discovery study with public summary statistics (v2.3) |
| `ms_independent_combined` | Generation Scotland ∪ HOLBAEK: supported if either supports; otherwise refuted if either refutes; otherwise indeterminate / not_assessed (v2.3) |
| `gs_multi_gene_entry` | True if the Generation Scotland entry maps to more than one gene (none among labelled CS in this release) (v2.3) |

A Chinese-language dictionary is provided as `FIELD_DICTIONARY_zh.md`.
