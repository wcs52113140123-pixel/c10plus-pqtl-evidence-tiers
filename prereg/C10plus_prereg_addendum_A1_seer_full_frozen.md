# C10+ 预注册增补 A1（Seer 2025 质谱全量汇总统计扩展）——冻结版
**冻结时间**：2026-11-03 12:25:37 UTC
**母预注册**：C10plus_preregistration_v2_frozen.md（SHA256 126f515f4734baef319629584521bc492e4a70730067f293c76964c2b62a89b0）。除本增补明确定义者外，一切定义、协变量、模型、阈值、停止规则沿用母预注册。
**随机种子**：20261103

## English summary
Addendum A1 extends the mass-spectrometry (MS) truth layer from the Suhre 2025 MSPA supplementary tables (169 evaluable credible sets) to the full Seer Proteograph GWAS summary statistics (GWAS Catalog GCST90570713–GCST90572692; 1,980 protein-group studies; Tarkin/MGB discovery cohort, N=1,252; GRCh37; inverse-normal-scaled traits). It defines H3i-ext, H3ii-ext and H4-ext with the same endpoints, thresholds and stopping rules as H3i/H3ii/H4, and a pre-specified MS "support / powered non-replication" rule for each UKB-PPP cis credible set: index variant = highest-PIP 95%-CS variant present in Seer (MAF≥0.01); supported = Seer P < 0.05/K (K = number of evaluated credible-set × study pairs) with UKB-concordant direction; refuted = ≥99% power (expected Seer z = |z_UKB|·√(N_Seer/N_UKB)) and not supported; otherwise indeterminate. Coordinates: UKB variants are lifted GRCh38→GRCh37 (UCSC hg38ToHg19 chain, pyliftover) and matched on position and unordered alleles; the 0- vs 1-based offset is fixed empirically from allele-matched overlap rates before any association value is inspected. **This addendum is written after the MSPA-based H3/H4 results were known (non-blind with respect to those results) but before any Seer full summary statistics were downloaded or viewed.** The original H3/H4 remain the pre-registered results and are reported unchanged; the extension is reported alongside, including a subset that excludes all credible sets already labelled by MSPA. The addendum family {H3i-ext, H3ii-ext, H4-ext} is Holm-adjusted at α=0.05.

## 0. 盲法状态
1. **已知**：母预注册中 H3i（OR 1.15，P=0.71，不成立）、H3ii（20% vs 52%，OR 0.23，P=0.092，部分成立）、H4（ΔAUC 0.10，P=0.135，部分成立）以及探索性结果（PAV 与 MSPA 支持负相关）均已知。因此本扩展**对这些结果非盲**；动机是 H3/H4 的质谱样本量不足（169 个可判定 cs）。
2. **未知**：写定本增补前，**未下载、未查看**任何 Seer 全量汇总统计的关联值。只查看了：GWAS Catalog v2 研究元数据（GCST90570713：性状 “GALC protein level (protein group normalized intensity)”、发现队列 MGBB N=1,252、replication 样本 325 不在该文件）、meta.yaml（GRCh37、标注为 0-based、HRC 填补、plink --linear、性状已对年龄/性别/10 个 PC 残差化后逆正态变换）和文件表头列名（chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, variant_id, n, protein_id1, protein_id2）。
3. **数据重叠**：MSPA（Suhre 2025 ST9/ST10）使用的 MS 数据来自同一研究（Tarkin + QMDiab），本扩展使用 Tarkin 发现队列全基因组汇总统计。因此扩展标签与 MSPA 标签**不独立**；为此预设“ext-new”子集（见 §4）。

## 1. 数据与抽取
- 范围：GCST90570713–GCST90572692 全部 1,980 个研究（每个 = 一个 MS 蛋白组 × nanoparticle）。
- 研究→基因：`protein_id1`（UniProt，取 “;”/“,” 分隔的第一个登录号去异构体后缀）→ Ensembl（OT gene.uniprot_swissprot，再 uniprot_all）；失败则用性状描述首个词作为 symbol 匹配；均失败记为未映射（报告计数）。
- 抽取窗口：同源基因 GRCh38 [start−1 Mb, end+1 Mb] 用 hg38ToHg19 链（pyliftover）分别转换两端；一端失败时以成功端 ± (基因长度 + 1 Mb) 补齐；两端均失败则不抽取（报告计数）。为对 0/1-based 不确定性留余量，窗口两端各再放宽 1 kb。
- 抽取在租用的云服务器（a rented cloud server）上流式完成（curl | gunzip | 窗口过滤 | gzip），只保留窗口内行；全文件不落盘。每个输出记录行数与 SHA256；同时校验 GWAS Catalog md5sum.txt 中的源文件 MD5（对流式读取的字节流计算 MD5）。

## 2. 坐标与等位基因协调
- UKB cs 变异（OT `qtl_cs_locus`，GRCh38）用 hg38ToHg19 链转换到 GRCh37；一对多或失败者剔除。
- 匹配键：(chrom, pos37, {A1,A2} 无序)。**0/1-based 偏移**：在查看任何 beta/P 之前，计算 offset ∈ {0, +1} 时 Seer 窗口变异与 OT 变异的等位基因匹配率，取匹配率高者（只用位置与等位基因，不用关联值）；该偏移对所有研究统一。
- 效应方向：Seer beta 对齐到 OT alt 等位基因（effect_allele==alt → +；==ref → −）。回文 SNP（A/T、C/G）且 Seer EAF∈[0.4, 0.6] 不用于方向判定（视为不可评估）。

## 3. 质谱支持判定规则（每个 UKB cis cs × 同源基因的每个 Seer 研究）
- 可评估变异：UKB 95% 集变异中在 Seer 匹配成功、Seer MAF ≥ 0.01、n ≥ 0.8 × 该研究最大 n 者。
- **指标变异**：可评估变异中 UKB PIP 最大者（先导可评估时即先导）。
- **Bonferroni 阈值**：α_ext = 0.05 / K，K = 全部可评估 (cs, 研究) 对的数目（对应 Suhre 的 0.05/392）。
- **期望 Seer z**：z_exp = |z_UKB| × √(N_Seer / N_UKB)。z_UKB 由 UKB 先导 −log10P 换算（cap 于 −log10P=300），N_UKB = OT 研究 nSamples（缺失用 sampleSize），N_Seer = 指标变异的 n。假设标准化效应相同（逆正态尺度），即“若 Olink 信号为真实丰度效应，MS 应见到的效应”。
- **功效**：power = Φ(z_exp − z_{1−α_ext/2}) + Φ(−z_exp − z_{1−α_ext/2})。
- **主标签 `ms_ext_primary`**：supported = 指标变异 Seer P < α_ext 且对齐后方向与 UKB β（该变异的 UKB locus β；先导用 qtl_cs.beta）一致；refuted = power ≥ 0.99 且未 supported（含显著但方向相反者，另行计数）；其余 = indeterminate（剔除）。
- 基因有多个 Seer 研究（多个 NP/蛋白组）：任一对 supported → supported；否则任一对 refuted → refuted；否则 indeterminate。
- 次标签：(S1) 名义标签：supported = P < 0.05 且同向；refuted = 在 α=0.05 下 power ≥ 0.99 且未 supported。(S2) 全可信集标签：supported = 任一可评估 CS 变异 P < α_ext / n_eval（该 cs 可评估变异数）且同向；refuted 同主标签功效规则。
- 无可评估变异或基因无 Seer 研究 → not_assessed。

## 4. 扩展终点（沿用母预注册定义、阈值与停止规则）
- **H3i-ext**：在 ms_ext_primary 可判定 cs 中，logistic `ms_ext ~ A_eqtl + mlp`（基因聚类 SE），预测 OR ≥ 1.5；报告 Fisher；阴性对照 W_gtex。
- **H3ii-ext**：concordant=1 中 CNP（A_eqtl=0 且 PAV=1）vs CA（A_eqtl=1）的 MS 支持率，Fisher 精确检验（双侧），预测 OR ≤ 0.5。停止规则：可判定 cs < 60 或 CNP/CA 任一组 < 10 → 不确定。
- **H4-ext**：目标 = ms_ext_primary；训练 = 奇数常染色体 + X，留出 = 偶数常染色体；特征、L2 logistic（C=1）、基线 B0/B1 与母预注册相同；主终点 DeLong 单侧 M_full > B0，阈值 ΔAUC ≥ 0.05。停止规则：留出可判定 < 40 或任一类 < 10 → 不确定。
- **多重比较**：增补族 {H3i-ext, H3ii-ext, H4-ext} Holm，α = 0.05；判定规则同母预注册（成立/部分成立/不成立/不确定）。
- **ext-new 子集（预设，报告但不另作确认性检验）**：剔除所有已有 MSPA 标签（ST10/ST9 匹配到的 cs，无论主标签是否可判定）后的 cs，重复 H3i-ext、H3ii-ext、H4-ext（H4-ext-new：在 ext-new 的留出集上评估同一模型）。若 ext 与 ext-new 方向相反，在论文中明确标出。
- **一致性 QC**：同时具有 MSPA 主标签与 ext 主标签的 cs，报告一致率与 Cohen κ。
- **次要（描述性）**：PAV、平台一致、A_iso 对 ms_ext 的 Fisher；Tier 1–4 的 ms_ext 支持率与趋势；S1、S2 标签下重复 H3i-ext/H3ii-ext；剔除 OT β 异常（flag_beta_inconsistent）cs 的敏感性。

## 5. 论文中的报告方式
- 母预注册 H3/H4 结果为原始确认性结果，原样报告，不被替换或改写。
- 扩展结果以“预注册增补（对 MSPA 结果非盲，对 Seer 全量数据盲）”单列，与原结果并排呈现（同一表的两列：MSPA-preregistered / Seer-full addendum，外加 ext-new）。
- 结论按两者共同权衡：一致则加强；不一致则同时报告并讨论（含 MS 数据重叠、功效、标签定义差异），不选择性强调。

## 6. 资源表新增字段
seer_ext_studies（GCST 列表）、seer_ext_index_variant、seer_ext_beta_aligned、seer_ext_se、seer_ext_p、seer_ext_power、seer_ext_zexp、ms_ext_primary、ms_ext_nominal、ms_ext_cswide、seer_ext_opposite_sig。

## 修订记录
（空）
