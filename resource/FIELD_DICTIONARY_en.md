# C10+ cis-pQTL evidence-tier resource v2.3: field dictionary

*English translation of `FIELD_DICTIONARY_zh.md` (v2.3). If the two versions differ, the Chinese version is authoritative.*

Files: `C10plus_cis_pqtl_robustness_grades_v2.3.parquet` / `.tsv.gz` (101 columns; the 76 v2.2 columns are unchanged value for value, and 25 addendum A4 columns were appended). 14,472 rows = all UKB-PPP European-ancestry cis credible sets (Open Targets 26.09, SuSiE-inf) in 1,879 genes. Definitions follow the pre-specified analysis plan v2.0 (SHA-256 126f515f…) and addenda A1, A2 and A4.

## Usage notes
- `tier` is a rule-based annotation from the pre-specified analysis plan; it does not depend on model training. **T4 means that affinity measurements alone cannot separate a true abundance effect from an epitope effect.** It is not a false-positive call (for example, functional missense variants in PCSK9 also fall into T4); such signals need mass-spectrometry or other orthogonal evidence.
- Tiers annotate PAV and transcript (eQTL) support. They are not a validated reliability ranking: in the analysis, transcript support did not predict mass-spectrometry confirmation.
- `platform_concordance` is informational only: cross-affinity-platform concordance was higher for PAV-containing signals, so for these signals concordance is not orthogonal validation. Overall, platform concordance was positively associated with MS support in the two independent MS cohorts (Generation Scotland, HOLBAEK; combined OR 2.92) but only weakly, and not significantly, in the two Seer layers.
- `p_ms_support_model` transfers poorly (held-out AUC 0.68; reversed when transferred to other truth labels). Do not use it alone for filtering.
- Colocalisation uses per-variant logBF within Open Targets 99% credible sets (not full-region lbf) and is an approximation.
- Credible sets on chromosome X have no MS support labels: none of the MS layers releases chrX summary statistics.

## Fields
| Field | Description |
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
| `dist_lead_to_tss` | \|lead − cognate TSS\| (bp) |
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
| `tier` | rule-based tier 1–4 (pre-specified analysis plan) |
| `tier_label` | T1 anchored no PAV; T2 anchored+PAV or isoform-only; T3 no anchor no PAV; T4 no anchor + PAV (regardless of platform concordance) |
| `p_ms_support_model` | L2 logistic (all chromosomes, MS target) predicted probability — low external validity, see H4 |
| `p_xa_rep_model` | L2 logistic (xa target) predicted probability |
| `model_confidence` | tertile of p_ms_support_model |
| `flag_concordant_unanchored_PAV` | platform-concordant but unanchored and PAV-containing (key caution class) |
| `lead_beta_to_expected_ratio` | \|beta\| / (z / sqrt(2pqN), N=34,000) |
| `flag_beta_inconsistent` | ratio>5 (OT beta inconsistent with P/MAF; 204 CS) — post hoc |

## Tier distribution
- Tier 1: 2539
- Tier 2: 682
- Tier 3: 10046
- Tier 4: 1205
- flag_concordant_unanchored_PAV = 339

## Fields added in v2 (addendum A1: Seer full summary statistics, GCST90570713–GCST90572692, Tarkin N≈1,252)
| Field | Description |
|---|---|
| `ms_ext_primary` | supported (index variant Seer P<0.05/1,779 with concordant direction) / refuted (power≥0.99 and not supported) / indeterminate / not_assessed |
| `ms_ext_nominal` | as above, threshold P<0.05, power computed at α=0.05 |
| `ms_ext_cswide` | supported if any evaluable 95% CS variant has P<α/n_eval with concordant direction |
| `seer_ext_opposite_sig` | a significant (cs, study) pair with opposite direction exists |
| `seer_ext_studies` | GCST studies evaluated |
| `seer_ext_index_variant` | index variant in the highest-powered study (the evaluable 95% CS variant with the largest PIP) |
| `seer_ext_beta_aligned` / `seer_ext_se` / `seer_ext_p` | Seer effect (aligned to the OT alt allele; inverse-normal scale), SE, P |
| `seer_ext_power` / `seer_ext_zexp` | power at the Bonferroni threshold; expected z = \|z_UKB\|·√(N_Seer/N_UKB) |

Note: the Seer files contain no chromosome X. GRCh37 coordinates were converted with hg38ToHg19 and matched on position plus unordered alleles (offset 0). MSPA and the Seer full statistics come from the same MS study; their primary labels agree in 94.6% (κ=0.89, n=148), so they are not independent validation.

## Fields added in v2.2 (addendum A2: ARIC, Zhang 2022; Fenland, Pietzner 2021; both SomaScan v4)
| Field | Description |
|---|---|
| `aric_ea_has_cis` / `aric_aa_has_cis` | cognate gene has any cis signal in ARIC EA/AA (ST3/ST6/ST10) |
| `rep_aric_ea` / `rep_aric_aa` | any ARIC signal variant (lead / conditionally independent / fine-mapped) lies in the UKB 95% set; NA without a cis signal |
| `rep_aric_ea_strict` | an ARIC-EA fine-mapped variant with PIP≥0.5 lies in the UKB 95% set |
| `dir_aric_ea` | effect direction when replicated (same/opposite/undetermined/NA) |
| `fenland_has_cis` / `rep_fenland` / `dir_fenland` | Fenland (name mapping, 93.0%) cis lead lies in the UKB 95% set; direction |

## Fields added in v2.3 (addendum A4: independent mass-spectrometry truth layers; Generation Scotland serum SWATH, N=15,035; HOLBAEK plasma Astral, N≈1,909)
| Field | Description |
|---|---|
| `ms_gs_primary` | supported (index variant GS P<0.05/1,103 with concordant direction) / refuted (power≥0.99 and not supported) / indeterminate / not_assessed |
| `ms_gs_nominal` / `ms_gs_cswide` | nominal-threshold label / credible-set-wide label (same rules as A1) |
| `gs_entries` | GS entries evaluated (UniProt groups; only entries that map to a single UKB cis gene) |
| `gs_index_variant` / `gs_beta_aligned` / `gs_se` / `gs_p` / `gs_power` / `gs_zexp` | index variant in the highest-powered entry, effect aligned to the OT alt allele, SE, P, power, expected z |
| `gs_opposite_sig` | a significant (cs, entry) pair with opposite direction exists |
| `gs_multi_gene_entry` | always False in the primary analysis (ambiguous multi-gene entries enter sensitivity analyses only) |
| `ms_hb_primary` / `ms_hb_nominal` / `ms_hb_cswide` | HOLBAEK labels (P<0.05/720); `not_evaluable` = no public summary statistics for this protein (only proteins with a pQTL were released, so selection bias exists) |
| `hb_entries` / `hb_index_variant` / `hb_beta_aligned` / `hb_se` / `hb_p` / `hb_power` / `hb_zexp` / `hb_opposite_sig` | as above (HOLBAEK, GRCh38 harmonised) |
| `hb_evaluable` | cognate gene has HOLBAEK discovery-set summary statistics |
| `ms_independent_combined` | GS ∪ HOLBAEK: any supported → supported; otherwise any refuted → refuted |

Note: GS covers autosomes only (F9, SERPINA7 and CFP are not evaluable). GS coordinates are GRCh37, matched after hg38ToHg19 conversion on position plus unordered alleles (offset 0); HOLBAEK is GRCh38 (offset 0). MS measures peptides and GS has no PAV-aware quantification, so GS support for a PAV-positive CS cannot by itself demonstrate an abundance effect. Primary-label agreement: GS vs Seer full statistics κ=−0.02 (n=73); GS vs HOLBAEK κ=0.57 (n=161).
