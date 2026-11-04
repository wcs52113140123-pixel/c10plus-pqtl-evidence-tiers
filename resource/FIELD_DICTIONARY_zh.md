# C10+ cis-pQTL 稳健性分级资源 v2.3 —— 字段说明

文件：`C10plus_cis_pqtl_robustness_grades_v2.3.parquet` / `.tsv.gz`；14,472 行 = UKB-PPP EUR 全部 cis 可信集（OT 26.09，SuSiE-inf），1,879 个基因。生成：2026-11-03；定义见预设分析计划 v2.0（SHA256 126f515f…）及增补 A1–A4。

## 使用建议
- `tier` 是预注册的规则分级，不依赖模型训练；**T4 表示仅凭亲和测量无法区分真实丰度效应与表位效应**，不是“假阳性”判定（例如 PCSK9 的功能性错义变异也会落入 T4），需要质谱或其他正交证据。
- Tier 注释的是 PAV 与转录（eQTL）支持，不是经过验证的可靠性排序：本研究中，转录支持并不能预测质谱确认。
- `platform_concordance` 只作为信息列：含 PAV 的信号跨亲和平台一致率偏高，所以对这类信号，一致不是正交验证。整体而言，平台一致在两个独立质谱队列（Generation Scotland、HOLBAEK）中与质谱支持正相关（合并 OR 2.92），在 Seer 两层中关联弱且不显著。
- `p_ms_support_model` 的外部可迁移性差（留出 AUC 0.68，迁移到其他真值标签反向），不建议单独用于筛选。
- 共定位基于 OT 99% 可信集内逐变异 logBF（非全区域 lbf），属近似。
- X 染色体上的可信集没有质谱支持标签：所有质谱层都没有发布 chrX 汇总统计。

## 字段
| 字段 | 说明 |
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

## 分级分布
- Tier 1: 2539
- Tier 2: 682
- Tier 3: 10046
- Tier 4: 1205
- flag_concordant_unanchored_PAV = 339

## v2 新增字段（预注册增补 A1：Seer 全量汇总统计，GCST90570713–GCST90572692，Tarkin N≈1,252）
| 字段 | 说明 |
|---|---|
| `ms_ext_primary` | supported（指标变异 Seer P<0.05/1,779 且方向一致）/ refuted（功效≥0.99 且未支持）/ indeterminate / not_assessed |
| `ms_ext_nominal` | 同上，阈值 P<0.05，功效按 α=0.05 计算 |
| `ms_ext_cswide` | 任一可评估 95% CS 变异 P<α/n_eval 且同向为 supported |
| `seer_ext_opposite_sig` | 存在显著但方向相反的 (cs, 研究) 对 |
| `seer_ext_studies` | 参与评估的 GCST 研究 |
| `seer_ext_index_variant` | 功效最高的研究中的指标变异（PIP 最大的可评估 95% CS 变异） |
| `seer_ext_beta_aligned` / `seer_ext_se` / `seer_ext_p` | Seer 效应（对齐 OT alt 等位基因，逆正态尺度）、SE、P |
| `seer_ext_power` / `seer_ext_zexp` | 在 Bonferroni 阈值下的功效；期望 z = \|z_UKB\|·√(N_Seer/N_UKB) |

注：Seer 文件不含 X 染色体；GRCh37 坐标经 hg38ToHg19 转换后按位置+无序等位基因匹配，偏移 0。MSPA 与 Seer 全量来自同一 MS 研究，二者主标签一致率 94.6%（κ=0.89，n=148），不是独立验证。

## v2.2 新增字段（增补 A2：ARIC Zhang 2022；Fenland Pietzner 2021；均为 SomaScan v4）
| 字段 | 说明 |
|---|---|
| `aric_ea_has_cis` / `aric_aa_has_cis` | 同源基因在 ARIC EA/AA 有任一 cis 信号（ST3/ST6/ST10） |
| `rep_aric_ea` / `rep_aric_aa` | 任一 ARIC 信号变异（先导/条件独立/精细定位）落在 UKB 95% 集；无 cis 信号时为 NA |
| `rep_aric_ea_strict` | ARIC-EA 精细定位变异 PIP≥0.5 落在 UKB 95% 集 |
| `dir_aric_ea` | 复现时的效应方向（same/opposite/undetermined/NA） |
| `fenland_has_cis` / `rep_fenland` / `dir_fenland` | Fenland（名称映射，93.0%）cis 先导落在 UKB 95% 集；方向 |


## v2.3 新增字段（增补 A4：Generation Scotland Richmond 2026，血清 Scanning SWATH，N=15,035；HOLBAEK Niu 2025，血浆 Orbitrap Astral DIA，N=1,909）
| 字段 | 说明 |
|---|---|
| `ms_gs_primary` / `ms_gs_nominal` / `ms_gs_cswide` | GS 质谱标签（规则同 A1：主标签 K=1,103，α=4.5×10⁻⁵；名义；CS 全域） |
| `gs_opposite_sig` / `gs_entries` / `gs_index_variant` | 显著反向对；评估的 DataShare 条目；指标变异 |
| `gs_beta_aligned` / `gs_se` / `gs_p` / `gs_power` / `gs_zexp` | GS 效应（对齐 OT alt 等位基因）、SE、P、功效、期望 z |
| `ms_hb_primary` / `ms_hb_nominal` / `ms_hb_cswide` | HOLBAEK 质谱标签（主标签 K=720，α=6.9×10⁻⁵）；无公开汇总统计的蛋白为 not_evaluable |
| `hb_opposite_sig` / `hb_entries` / `hb_index_variant` / `hb_beta_aligned` / `hb_se` / `hb_p` / `hb_power` / `hb_zexp` | 同上，HOLBAEK |
| `hb_evaluable` | 该蛋白有 HOLBAEK 发现集公开汇总统计 |
| `ms_independent_combined` | GS ∪ HOLBAEK：任一支持即 supported；否则任一否定即 refuted |
| `gs_multi_gene_entry` | GS 条目映射到多个基因（本版带标签 CS 中无） |

注：GS 与 HOLBAEK 彼此一致（κ=0.57），但与 Seer 标签几乎不一致（κ −0.02 至 0.11）；四个 MS 层均不含 X 染色体。Tier 是证据模式注释，不是可靠性排序。
