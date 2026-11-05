# C10+ 预注册增补 A4：独立质谱真值层（Generation Scotland 为主，HOLBAEK 为次）
版本：A4 v1.0（冻结时间与 SHA256 见同名 .sha256 文件）
母预注册：v2.0（SHA256 126f515f4734baef319629584521bc492e4a70730067f293c76964c2b62a89b0）；增补 A1（SHA256 6d54bd54…035e2）、A2、A3。除本增补明确定义者外，一切定义、协变量、模型、阈值、判定规则沿用母预注册与 A1。
随机种子：20261103。资源表基线：v2.2（studyLocusId 键）；分级规则不改。

## English summary
Addendum A4 adds two mass-spectrometry (MS) pQTL truth layers that are sample- and technology-independent of Suhre 2025 (Seer/MSPA): Generation Scotland (GS; serum Scanning SWATH; N=15,035; primary) and HOLBAEK (Niu 2025; plasma Orbitrap Astral DIA; N≈1,909; secondary). MS support labels are computed per UKB-PPP cis 95% credible set with exactly the A1 rule (index = highest-PIP evaluable CS variant; supported = P<0.05/K with UKB-concordant direction; refuted = power ≥0.99 under equal standardized effect and not supported; else indeterminate), with power recomputed for each layer's N. Confirmatory family (Holm α=0.05): H3i-A4, H3ii-A4, H4-A4 on GS. H4-A4 applies the frozen A1 model (trained on Seer labels, odd chromosomes; not retrained) to GS labels on held-out even chromosomes. HOLBAEK endpoints are a secondary replication; HOLBAEK proteins without public summary statistics are "not evaluable" (selection bias declared). The author is non-blind to H3/A1 results.

## 0. 盲态声明
1. **对既往结果非盲**：作者已知 H3i（OR 1.15，P=0.71）、H3ii（20% vs 52%，P=0.092）、H4（AUC 0.68 vs 0.58，P=0.135）、探索性 PAV→MSPA 支持（34% vs 56%，P=0.012），以及 A1 的 H3i-ext（OR 1.21，P=0.50）、H3ii-ext（45% vs 46%）、H4-ext（AUC 0.636 vs 0.469，Holm 0.007）、κ(MSPA, Seer)=0.89。
2. **可行性阶段已接触的数据（如实列出）**：
   - GS：读取 DataShare README 与 GWAS_sumstats.zip 的 ZIP 中心目录（条目名、大小、偏移）。为确认格式，解压了 1 个条目（P02649_E7ERP7_H0Y7L5_E9PEV4，APOE 蛋白组）开头约 400 KB，只看表头与 chr10 起始若干行，**未查看任何 cis 区域**。
   - HOLBAEK：下载补充表（含 ST4 主 pQTL 列表与肽段证据分层）与 1 个 sumstats 的 yaml；正文仅做吞吐测速（约 60 MB 字节流，未解析）。
   - G&H Seer、GNHS：下载补充表，仅用于计算基因覆盖，**未与 UKB 比对**。G&H 与 GNHS 不进入本增补的任何终点。
   - 以上均未与 UKB-PPP CS 做变异层面匹配，未计算任何复现、方向或支持结局。
3. 覆盖估计（可行性报告）：GS 99 个 UKB cis 基因/1,093 CS；HOLBAEK 178 个基因/2,156 CS。这些是上限，不是可判定数。

## 1. 数据
- **GS（主层）**：Richmond 等，Nat Commun 2026，doi:10.1038/s41467-026-76474-8，PMID 42860425；Edinburgh DataShare doi:10.7488/ds/8011（hdl 10283/9082），CC BY 4.0；GWAS_sumstats.zip 内 439 个 `<UniProt组>.fastGWA.gz`；GRCh37（HRC）；列 CHR SNP POS A1 A2 N AF1 BETA SE P；fastGWA 约定 A1 为效应等位基因、AF1 为 A1 频率（抽取前以 README 核对；若 README 另有说明，以 README 为准并记录，在查看任何关联值之前完成）。
- **HOLBAEK（次层）**：Niu 等，Nat Genet 2025，doi:10.1038/s41588-025-02089-2，PMID 39972214；GWAS Catalog GCST90452968–GCST90453410（发现集 443 个研究，仅限有 pQTL 的蛋白）；使用 harmonised GRCh38 文件（含 .tbi）。复制集（GCST90453411 起）不使用。
- **选择偏倚声明**：HOLBAEK 只公开至少有一个 pQTL 的蛋白；测到但未公开的蛋白归为 `not_evaluable`，不能判为“未复现”。由于公开蛋白都有某种 pQTL，HOLBAEK 支持率可能被抬高，且两组比较（H3ii）在公开蛋白内部进行。

## 2. 条目→基因映射
- GS：条目名按“_”拆分为 UniProt 登录号（去异构体后缀），经 OT 26.09 gene.uniprot_swissprot，再 uniprot_all 映射到 Ensembl 基因；条目的基因集 = 所有登录号映射基因 ∩ UKB-PPP cis 基因（资源 v2.2）。
  - 主分析：恰好映射到 1 个 UKB cis 基因的条目（可行性估计 206 个）。
  - 映射到 ≥2 个 UKB cis 基因的条目（估计 4 个，如 P08519_P00747 = LPA/PLG）为歧义蛋白组，主分析排除；敏感性分析中对每个基因各计一次。
- HOLBAEK：研究→蛋白以 GWAS Catalog 研究性状/UniProt 注释为准，再按同一 OT 规则映射；映射失败计数报告。
- 同一基因有多个条目/研究：按 A1 的合并规则（任一 supported → supported；否则任一 refuted → refuted；否则 indeterminate）。

## 3. 抽取（冻结后执行）
- 窗口：同源基因 GRCh38 [start−1 Mb, end+1 Mb]。GS 用 A1 的 hg38ToHg19 转换规则（pyliftover，一端失败时以成功端 ±（基因长度+1 Mb）补齐），两端各放宽 1 kb；HOLBAEK 直接用 GRCh38 窗口。
- GS：以 HTTP Range 按 ZIP 本地文件头偏移只取所需条目字节（约 210 个条目、约 34 GB），流式解压后只保留窗口内行；全文件不落盘。校验：每个条目流式字节的 CRC-32 与 ZIP 中心目录记录比对（不一致则重取，二次失败记为抽取失败）；记录条目字节 SHA256、输出行数与输出文件 SHA256。若 DataShare 提供整包 MD5，因只取部分字节而无法核对整包，如实声明。
- HOLBAEK：优先 tabix 远程按窗口查询；随机抽 5 个研究整文件下载，核对 GWAS Catalog md5，并要求整文件按窗口过滤的结果与远程 tabix 结果逐行一致，否则改为整文件下载后过滤。
- 执行位置：本地 Linux 工作站（a local Linux workstation）；仅在明显更快时用租用的云服务器（a rented cloud server；用后清理工作目录、不开放端口）；不使用其他服务器。

## 4. 坐标与等位基因协调（同 A1）
- UKB CS 变异（OT qtl_cs_locus，95% CS）：GS 用 GRCh37 坐标（hg38ToHg19 转换，一对多或失败者剔除）；HOLBAEK 用 GRCh38 坐标。
- 匹配键：(chrom, pos, {A1,A2} 无序)。0/1-based 偏移：在查看任何 beta/P 之前，分别对 GS 与 HOLBAEK 计算 offset ∈ {0, +1} 的等位基因匹配率，取高者，全层统一。
- 方向：效应对齐到 OT alt 等位基因；回文 SNP 且 EAF∈[0.4,0.6] 不用于判定。

## 5. 质谱支持判定（每个 UKB cis CS × 同源基因的每个条目/研究；完全沿用 A1 §3，按层分别计算）
- 可评估变异：95% CS 变异中匹配成功、MAF≥0.01、n≥0.8×该条目最大 n（HOLBAEK 无 n 列时 n 取 1,909，即研究元数据样本量）。
- 指标变异 = 可评估变异中 UKB PIP 最大者。
- α_layer = 0.05/K_layer（K = 该层全部可评估 (CS, 条目) 对数）。
- z_exp = |z_UKB|·√(N_layer/N_UKB)；power = Φ(z_exp − z_{1−α/2}) + Φ(−z_exp − z_{1−α/2})。
- supported = 指标变异 P<α 且方向与 UKB 一致（P 记为 0 者视为 <α）；refuted = power≥0.99 且未 supported（显著但反向者另计）；否则 indeterminate。
- 次标签：名义标签（α=0.05）、全可信集标签（同 A1 S1/S2）。
- 无可评估变异或基因无条目/研究 → not_assessed（HOLBAEK 无公开汇总统计 → not_evaluable）。

## 6. 终点
### 6.1 确认性族（GS，Holm α=0.05）
- **H3i-A4**：GS 主标签可判定 CS 中，logistic `ms_gs ~ A_eqtl + mlp`（基因聚类 SE），预测 OR≥1.5；同时报告 Fisher 与阴性对照 W_gtex。
- **H3ii-A4**：concordant=1 中 CNP（A_eqtl=0 且 PAV=1）vs CA（A_eqtl=1）的 GS 支持率，Fisher 双侧，预测 OR≤0.5。停止：可判定 CS<60 或任一组<10 → 不确定。
- **H4-A4**：使用**冻结的 A1 模型**（FULL、B0、B1，各自在 A1 训练集——奇数常染色体+X、目标 ms_ext_primary——上拟合；不重训。因 L2 logistic 确定性，按 A1 代码在相同训练数据上复现拟合，并核对标准化系数与 A1 记录一致，误差 <1e-6）。在偶数常染色体上 GS 主标签可判定的 CS 上评估；主终点 DeLong 单侧 FULL > B0，阈值 ΔAUC≥0.05 且 Holm P<0.05；同时报告 vs B1（平台一致+功效）。停止：评估集<40 或任一类<10 → 不确定。次要：全染色体 GS 可判定 CS（奇数染色体与 A1 训练特征重叠，标签不同，仅描述）。
- **全层停止规则**：GS 主标签可判定 CS<100 → GS 三个确认性终点均判“不确定”。
- 判定规则同母预注册（成立 / 部分成立 / 不成立 / 不确定）。

### 6.2 HOLBAEK 次层复制（非确认性；族内 Holm 仅作参考）
- 以 HOLBAEK 主标签重复 H3i、H3ii、H4（同上定义、阈值与停止规则，停止规则的总数门槛同为 100）。

### 6.3 次要终点（均预先声明；非确认性）
- **PAV→MS 支持**（原为探索性，此处作为独立复现的次要终点）：GS 与 HOLBAEK 分别做 PAV（cognate PAV in 95% CS，VEP116 all-transcript）vs 非 PAV 的支持率 Fisher 双侧；预设方向 OR<1。
- **与 Seer 标签一致性**：GS vs ms_ext_primary（A1）、GS vs ms_primary（MSPA）、GS vs HOLBAEK、HOLBAEK vs ms_ext_primary 的一致率与 Cohen κ。
- **GS 新基因子集**：A1 中无 Seer 研究（seer_ext_studies 为空）的基因，重复 H3i/H3ii（描述性）。
- Tier 1–4 的 GS/HOLBAEK 支持率与 Cochran–Armitage 趋势；flag_concordant_unanchored_PAV CS 的支持率。
- 名义标签、全可信集标签下重复 H3i/H3ii；剔除 flag_beta_inconsistent 的敏感性。
- 歧义多基因条目纳入的敏感性。
- 合并标签 GS∪HOLBAEK（任一 supported → supported；否则任一 refuted → refuted）重复 H3i/H3ii（描述性）。

## 7. PAV 阳性 CS 中 MS 支持的解释规则
MS 测的是肽段。如果 PAV 改变了某个被用于定量的肽段，携带变体等位基因者的该肽段不被参考谱库识别，蛋白强度下降。这会在 MS 侧产生与亲和平台“表位效应”同向的伪影：两平台都显示变体等位基因降低信号，从而出现“假支持”。GS 没有 PAV 感知定量（未剔除含变体肽段），因此：
- 主分析不改标签，但在解释中注明：PAV 阳性 CS 的 GS 支持不能单独证明丰度效应。
- **敏感性 S-PAVdir（预设）**：PAV 阳性 CS 中，取 CS 内 PIP 最高的 cognate PAV；若该 PAV alt 等位基因的 UKB β<0（Olink 下降），则 GS“supported”在敏感性中改为 indeterminate（无法排除肽段丢失伪影）；若 UKB β>0（升高），肽段丢失不能解释同向升高，保留 supported。重复 H3ii-A4 与 PAV→MS 支持。
- **HOLBAEK 肽段分层敏感性**：若 ST4 的肽段证据分层可将 CS 对应的 HOLBAEK cis pQTL（同基因，ST4 lead 位于该 CS 或与 UKB 指标变异相同）标为可能受变体肽段影响的层级（按 Niu 2025 定义），则在敏感性中把这些 CS 的 HOLBAEK 标签设为 indeterminate，重复 PAV→MS 支持与 H3ii-HB。

## 8. 资源表 v2.3 新增字段
gs_entries、gs_index_variant、gs_beta_aligned、gs_se、gs_p、gs_power、gs_zexp、ms_gs_primary、ms_gs_nominal、ms_gs_cswide、gs_opposite_sig、gs_multi_gene_entry；hb_studies、hb_index_variant、hb_beta_aligned、hb_se、hb_p、hb_power、hb_zexp、ms_hb_primary、ms_hb_nominal、ms_hb_cswide、hb_opposite_sig、hb_evaluable；ms_independent_combined。

## 9. 报告方式
原 H3/H4 与 A1 结果原样保留；A4 与其并排报告（MSPA / Seer 全量 / GS / HOLBAEK 四列）。一致则加强，不一致则同时报告并讨论（血清 vs 血浆、技术、功效、PAV 肽段伪影、选择性公开），不选择性强调。

## 修订记录
（空）
