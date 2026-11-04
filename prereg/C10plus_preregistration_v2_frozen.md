# C10+ 全量分析预注册（v2.0，冻结版）
**题目（暂定）**：Cross-platform concordance is not orthogonal validation: transcript-anchored and mass-spectrometry-calibrated grading of plasma cis-pQTLs
**冻结时间**：2026-11-03 11:02:32 UTC（本文件写定后计算 SHA256；其后任何改动只能以“修订记录”追加，不得修改终点、阈值、检验或分组定义）
**数据**：thinkstation `/home/weichangsheng/.claude-science-scratch/c10_pqtl/v1`（DuckDB/parquet 23 张表；README v1，2026-11-03）
**随机种子**：20261103（所有抽样、匹配、交叉验证、置换）

---
## English summary (abstract of the plan)
We test whether agreement between two affinity proteomics platforms (Olink in UKB-PPP vs SomaScan in INTERVAL/Sun 2018 and deCODE/Eldjarn 2023) is a valid orthogonal validation of cis-pQTLs, or whether transcript-level anchoring (SuSiE-based colocalisation with same-gene eQTL) and antibody-independent mass spectrometry (MS) are more reliable evidence. Units are UKB-PPP EUR cis credible sets (Open Targets 26.09, SuSiE-inf; 14,472 credible sets, 1,879 genes). Five hypotheses: **H1** (internal reproduction, Sun 2018; overlaps the pilot data and is NOT confirmatory) – eQTL anchoring predicts cross-affinity replication beyond a wrong-gene negative control and power matching; **H2** (confirmatory, deCODE SomaScan v4 lead/LD-class signals in 888 genes not assayed by Sun 2018) – the same anchor effect, plus a power-matched test that protein-altering variants (PAV) increase binary cross-affinity replication, with effect-correlation/direction analyses to reconcile with Eldjarn 2023; **H3** (confirmatory, core) – in the MS truth layer (Suhre 2025 MSPA ST10/ST9; Seer MS pQTL ST1–ST3 as secondary), (i) anchored signals are MS-supported more often, and (ii) affinity-concordant but unanchored PAV-containing signals are MS-supported less often than concordant anchored signals; **H4** – a robustness grade trained on odd chromosomes predicts MS support and cross-ancestry replication on held-out even chromosomes better than platform concordance alone; **H5** – anchored signals replicate more often (same direction, P<0.05) in pooled non-European UKB-PPP ancestries (AFR/CSA/EAS/MID/AMR). Confirmatory primary endpoints (H2a, H2b, H3i, H3ii, H4, H5) are Holm-adjusted at family-wise α=0.05; each requires both an adjusted P<0.05 and a pre-specified effect-size threshold. Gene-clustered robust SEs are used throughout; stopping/downgrade rules declare under-powered tests inconclusive rather than negative.

---
## 0. 盲法状态与已知信息披露（写定前已看到的东西）
1. **预实验（2026-11-02）**已在 Sun 2018 的 424 个基因（4,153 个 UKB cis 可信集）上看到 G1–G3 和 PAV 结果（锚定 OR 3.78；PAV OR 2.26；四格复现率 7.0/20.0/37.2/50.0%）。因此 **H1 只是预实验的全量重做/内部复现，不是确认性证据**；H2b、H3(ii) 的方向预测由该预实验产生。
2. 写定前查看过：各表结构（列名、行数）；ext_pqtl_layers 各层边际计数；覆盖矩阵（v3/v4/MS 检测的基因数）；Suhre 2025 MSPA 补充表 ST9/ST10/ST12/ST13 的列名与真值列边际计数（如 ST10：power>99% 167/392，MSPA>0.8 64/392）。**未计算任何锚定/PAV/平台一致 × 质谱/deCODE/非欧祖先结局的交叉表**，也未对 deCODE、质谱、非欧数据做任何 QTL×QTL 共定位。
3. 已知的已发表结论：Eldjarn 2023（PAV 降低平台间效应相关；eQTL/PAV 的 LD 分类）；Suhre 2025（两亲和平台同时检出被视为无表位效应的强指标；约 30% 亲和 pQTL 被 MS 支持）；Nicholas 2026（平台一致的 cis-pQTL 更常含 GTEx eQTL 先导变异）。

## 1. 通用定义
- **分析单位**：OT 26.09 `qtl_cs` 中 `projectId='UKB_PPP_EUR'` 且 `isTransQtl=false` 的 95% 可信集（cs）。cs 的变异集合 = `qtl_cs_locus` 中该 cs 的全部变异（OT 只给 99% 集；`is95` 标 95% 集）。
- **QTL×QTL 共定位**：沿用预实验方法——用两侧 SuSiE 逐变异 logBF 按 coloc 单效应公式（Giambartolomei 2014；Wallace 2021 coloc.bf_bf 思路）计算 PP.H4；一侧缺失的变异 logBF 记为 0。先验 p1=p2=1e-4，p12=1e-5。**共定位阳性 = PP.H4 ≥ 0.8**（敏感性 0.9/0.95）。对一个 UKB cs 与一组对侧 cs，取最大 H4。
- **锚定（anchor）**：
  - `A_gtex`：与同基因 OT GTEx-v10 eQTL cs（任一组织）共定位。
  - `A_nongtex`：与同基因 OT eQTL Catalogue eQTL 或 sc-eQTL cs（projectId≠GTEx-v10）共定位。
  - `A_eqtl` = A_gtex ∨ A_nongtex（H3–H5 的主暴露）。
  - `A_isoform`：与同基因任一 sQTL 或 tuQTL cs（GTEx-v10 或 eQTL Catalogue）共定位。
  - 注意：eQTL Catalogue 的 `aptamer` 数据集（Sun 2018 SomaScan）属复现层，**不作锚定**；OT 的 Sun_2018 是 pQTL，不在 eqtl 中。
- **错基因阴性对照** `W_gtex`：与任一**非同源基因**的 GTEx-v10 eQTL cs 共定位（候选 = 与该 UKB cs 至少共享 1 个变异的非同源基因 eQTL cs）。
- **PAV**（主定义）：`pqtl_cs_pav_alltx_vep116.cs95_pav_cognate_alltx_vep116`（VEP 116 全转录本，同源基因，95% 集内含 missense/stop_gained/stop_lost/start_lost/frameshift/inframe_ins/del/protein_altering/splice_donor/splice_acceptor）。敏感性：OT 规范转录本 `cs95_pav_cognate_canonical`；PIP 合计 ≥0.2（预实验定义 `pp_sum_pav_cognate_alltx_vep116≥0.2`）；先导变异本身为 PAV（`pqtl_cs_leadpav`）。
- **功效协变量**：`mlp` = min(UKB lead −log10P, 300)/10；`logn95` = log(max(95% 集大小,1))；`absb` = |lead β|；`maf` = min(eaf, 1−eaf)。
- **基因聚类推断**：所有 logistic 回归报告按基因聚类的稳健 SE（statsmodels GLM cov_type='cluster'）作为主推断；同时报告模型 SE。
- **HLA 区**：chr6:28,510,120–33,480,577（GRCh38）。

## 2. 复现层定义
### 2.1 Sun 2018（SomaScan v3，INTERVAL，N≈3,301；H1）
`rep_sun` = UKB cs 与同基因 OT Sun_2018 cis cs 共定位（H4≥0.8）。分母：Sun 2018 检测了该基因（`assayed_SomaScan_v3_Sun2018`，424 基因）。协变量 `sun_has_cis` = 该基因在 Sun 2018 有任一 cis cs。

### 2.2 deCODE SomaScan v4（Eldjarn 2023 ST19/ST33；H2）
- 只有先导信号与 LD 类（LD class，Eldjarn 定义为与先导 r² 高的变异），无 PIP。**信号匹配规则（主）**：同基因的 deCODE cis 信号（ST19 `cis_trans_src='cis'` 或 ST33 cis，或 `is_cis_1mb`）中，其先导变异或 LD 类成员（`cs_members_src` 中的 chr:pos 直接用 GRCh38 位置；rsID 经 `variant_core.rsIds` 映射为位置）与 UKB cs 的 **95% 集变异按 (chrom,pos) 有任一重合** → `rep_decode=1`。
- 敏感性：(a) 只用 deCODE 先导变异落入 UKB 99% 集；(b) 允许 UKB 99% 集。
- 分母（H2 基因集）：`assayed_SomaScan_v4_deCODE=True` 且 `assayed_SomaScan_v3_Sun2018=False` 的 UKB cis 基因（888 个）。注意：v4“是否检测”由 Eldjarn ST19/ST33 中是否出现该基因（cis 或 trans 信号）推断，属近似。协变量 `decode_has_cis`。

### 2.3 平台一致（concordance，用于 H3/H4）
`concordant` = rep_sun ∨ rep_decode；`discordant` = 该基因在 v3 或 v4 被检测但两者都未复现；`na` = 基因未在任何 SomaScan 层检测。敏感性：一致须同时满足效应方向相同（可判定时）。

## 3. 真值层定义
### 3.1 质谱 MSPA（Suhre 2025 Nat Genet，doi:10.1038/s41588-025-02413-w，PMID 41310232；ST10 = 392 个 UKB-PPP pQTL，ST9 = 319 个 deCODE pQTL）
- 匹配：MSPA 行的 UniProt → Ensembl 基因（gene.uniprot_swissprot/uniprot_all），变异 hg38 位置（ST10 `GENPOS (hg38)`；ST9 `SNP (build 38)`）落入同基因 UKB cis cs 的 99% 集（按 chrom:pos；多 cs 命中取该变异 PP 最大的 cs）。一个 cs 有多行时优先 ST10，再取 power 最高的行。未匹配行计数并报告。
- **主标签** `ms_primary`：supported = Bonferroni 复现（ST10 `replicated (p<0.05/392)`；ST9 `replicated (p<0.05/319)`）；refuted = `power>99%` 且未复现；其余（功效不足且未复现）= 不确定，剔除。
- 次标签 `ms_mspa`：supported = MSPA>0.8；refuted = MSPA<0.2；中间值剔除。
- 三级标签 `ms_seer`（Seer MS pQTL 发现层，ST1/ST2/ST3）：对 MS 检测了该蛋白（ST1 出现）的基因上的 UKB cs，若任一 Seer cis 信号（ST2 先导、ST3 条件信号，或 ST1 每蛋白最显著变异且 −log10P ≥ 7.3）的变异按 allele_key 落入 UKB 99% 集 → supported，否则 not supported；分析时调整 UKB 功效（mlp、absb、maf）。
- 已知局限：MS 肽段若覆盖 PAV 氨基酸也会产生肽段假象；Suhre 的 MSPA 已在肽段层面处理（PEPanno、PAV rsID 列），我们另做敏感性：剔除 `MBRplus.PEPanno.PAVrsid` 非空的行。

### 3.2 非欧祖先（UKB-PPP Synapse syn51364943，AFR/CSA/EAS/MID/AMR，cis ±1 Mb 全汇总统计；H4/H5）
- 指标变异：EUR 先导变异；若某祖先中缺失、MAF<0.01 或 INFO<0.8，则改用该 cs 中通过过滤、EUR PP 最大的变异；均不通过则该祖先缺失。
- 等位基因协调：OT `variant_id` 的 alt 为效应等位基因；非欧 `allele1` 为效应等位基因（regenie），allele1==alt → β 不变，allele1==ref → 取反，否则剔除。EUR β 用该变异在 UKB locus 的 β（先导用 qtl_cs.beta）。
- 跨祖先合并：对可用祖先做固定效应逆方差合并。
- **主结局** `xa_rep` = 合并 β 与 EUR β 同号 且 合并 P<0.05。次要：同号率；效应一致性斜率（合并 β 对 EUR β 的过原点回归）；分祖先（AFR、CSA 样本最大）。

## 4. 假设、终点与判定
**判定规则（全部假设通用）**：成立 = Holm 校正后 P<0.05（确认性主终点族）且效应量达阈值，且对应阴性对照不显著/显著更小；部分成立 = 仅满足其一，或主终点成立而关键次要/阴性对照不支持；不成立 = 都不满足或方向相反；不确定 = 触发停止/降级规则（功效不足）。
**确认性主终点族（Holm，α=0.05）**：H2a、H2b、H3i、H3ii、H4、H5（6 项）。H1 不在确认族内（内部复现，用预实验门槛独立判定）。次要终点在各假设内 Holm。

### H1（内部复现，Sun 2018；非确认性；数据与预实验重叠）
- 模型：logistic `rep_sun ~ A_gtex + mlp + logn95 + PAV + sun_has_cis`，基因聚类 SE。
- 门槛（同预实验）：G1 OR≥1.5 且 P<0.001，且 `W_gtex` 的 OR 显著小于 `A_gtex`（联合模型差值单侧 P<0.05）；G2 `A_nongtex` OR≥1.5 且 P<0.01；G3 功效匹配子集（UKB P<1e-30 且 sun_has_cis=1）OR≥1.3 且 P<0.05。
- 次要：PAV 效应（调整 A_gtex 与功效）；PAV 功效匹配（同 H2b 方法）；四格复现率；效应相关与方向（同 H2c）。

### H2（确认性，deCODE v4，888 个 Sun 2018 未覆盖基因）
- **H2a 主终点**：logistic `rep_decode ~ A_gtex + mlp + logn95 + PAV + decode_has_cis`，基因聚类 SE；预测 OR≥1.5（阈值），单侧方向 >1；检验为双侧 Wald P。阴性对照：`W_gtex` 在同模型 OR 显著小于 `A_gtex`（联合模型差值单侧 P<0.05）。次要：`A_nongtex`（OR≥1.5）；功效匹配子集（UKB P<1e-30 且 decode_has_cis=1）OR≥1.3。
- **H2b 主终点（回应“PAV 复现高只因功效高”）**：PAV 功效匹配——每个 PAV+ cs 在 PAV− cs 中按 (mlp, absb, maf, logn95) 的 Mahalanobis 距离做 1:3 最近邻匹配（无放回，卡钳 0.5 SD 于马氏距离标准化后，限定同一 decode_has_cis 层），在匹配样本上 logistic `rep_decode ~ PAV`，按匹配集聚类 SE。预测 OR>1，阈值 OR≥1.3。同时报告：全样本协变量调整（四个功效变量用 4 df 自然样条）的 PAV OR；匹配前后标准化均值差（SMD<0.1 视为平衡）。
- **H2c 次要（调和 Eldjarn 2023）**：在 deCODE 先导变异位于 UKB locus 且等位基因可协调的 cs 上，比较 PAV+ 与 PAV− 的 (i) 二元复现率，(ii) 平台间 β 的 Pearson 相关（Fisher z 差值检验，预测 PAV+ 更低），(iii) 方向相反比例（Fisher 精确检验，预测 PAV+ 更高），(iv) 只在同向对中计算的相关。H1 对 Sun 2018 做同样分析（共享变异取 PP_ukb×PP_sun 最大者）。**预测**：PAV 同时提高二元复现、降低效应相关/提高反向比例——二者不矛盾（衡量的量不同）。

### H3（确认性，核心：质谱真值）
- **H3i 主终点**：在 `ms_primary` 可判定的 cs 中，MS 支持率 A_eqtl=1 vs 0；主检验 logistic `ms_supported ~ A_eqtl + mlp`（基因聚类 SE），预测 OR≥1.5；同时报告 Fisher 精确检验与原始比例。阴性对照：`W_gtex` 替换 A_eqtl，预期 OR 不显著。
- **H3ii 主终点（“平台一致≠正交验证”）**：在 concordant=1 的 cs 中比较 组 CNP（无锚定 A_eqtl=0 且 PAV=1）与 组 CA（A_eqtl=1，PAV 任意）的 MS 支持率；主检验 Fisher 精确检验（双侧），预测 OR(CNP vs CA)≤0.5（阈值）。次要：Firth/普通 logistic 调整 mlp；同样比较 CN0（无锚定无 PAV）vs CA；以及“discordant 且有锚定” vs CNP。
- 次要：`ms_mspa`、`ms_seer` 两个次标签下重复 H3i/H3ii；剔除 PAV 肽段行的敏感性。
- **停止/降级规则**：若 ms_primary 可判定 cs <60，或 H3ii 的 CNP 或 CA 组 <10 个可判定 cs，则 H3ii 主终点判为“不确定（功效不足）”，其余标签只作描述；不得以次标签替代主终点。

### H4（稳健性分级的建模与留出验证）
- **划分**：训练 = 奇数常染色体 + X；留出 = 偶数常染色体（固定，不调换）。
- **特征（全模型 M_full）**：A_eqtl、log1p(锚定数据集数)、A_isoform、PAV、pp_sum_pav、PAV×(1−A_eqtl)、平台一致三分类（concordant/discordant/na 哑变量）、mlp、absb、maf、logn95、purityMinR2、log10(1+|lead 到 TSS 距离|)。标准化后 L2 logistic（C=1.0，不调参）。
- **基线**：B0 = 仅平台一致三分类；B1 = 平台一致 + mlp + absb + maf。
- **目标与主终点**：(a) **H4 主终点**：以 `ms_primary` 为目标，在训练染色体上拟合 M_full 与 B0，在留出染色体上比较 AUC，DeLong 单侧检验 M_full>B0，阈值 ΔAUC≥0.05。(b) 次要：同法对 `xa_rep`（H5 结局）训练/留出评估；MS 训练模型对留出 xa_rep 的迁移 AUC；M_full vs B1；校准（留出集校准截距、斜率、Brier）；标签置换零分布（1,000 次）。
- **规则分级（预设，不训练）**：
  - Tier 1 高：A_eqtl=1 且 PAV=0；
  - Tier 2 中：A_eqtl=1 且 PAV=1，或 A_eqtl=0、A_isoform=1、PAV=0；
  - Tier 3 未定：A_eqtl=0、A_isoform=0、PAV=0（含“组织/表达未覆盖”可能）；
  - Tier 4 警示：A_eqtl=0 且 PAV=1（含 A_isoform=1 且 PAV=1）——**无论平台是否一致**。
  - MS 状态（supported/refuted）与跨祖先（复现/不复现/不可评估）作为独立标注列，不改 Tier；置信度 = 全染色体重拟合 M_full（MS 目标）预测概率。
  - 次要：Tier 1→4 的 MS 支持率与 xa_rep 趋势（Cochran–Armitage）。
- **停止规则**：留出集 ms_primary 可判定 <40 或任一类 <10 → H4 主终点“不确定”，以 xa_rep 结果作描述。

### H5（跨祖先）
- **主终点**：logistic `xa_rep ~ A_eqtl + mlp + absb + maf_eur + log(se_meta) + logn95`，基因聚类 SE；预测 OR≥1.3。`log(se_meta)` 吸收非欧样本量与等位基因频率造成的功效差。
- 次要：效应一致性斜率按 A_eqtl 分层（交互项检验）；分祖先（AFR、CSA）；同号率；`W_gtex` 阴性对照。
- PAV 敏感性：按 PAV 分层；剔除 |MAF_nonEUR − MAF_EUR| > 0.1 的 cs（等位频率差异）；只用 EUR 先导变异（不替换）。
- **停止规则**：可评估 cs <500 → 不确定。

## 5. 敏感性分析（预设，全部在 H1 G1、H2a、H2b、H3i、H5 主估计上重做）
1. H4 阈值 0.9、0.95；
2. coloc 先验 p12 = 5e-6、1e-6；以及只用两侧共有变异的 coloc（bf_bf 交集版）；
3. 位点窗口：只保留先导距 TSS ≤1 Mb、≤500 kb 的 cs；
4. 每基因只取最强 cs（lead only）；
5. 剔除 HLA 区；
6. 剔除 trans 混入：剔除其 95% 集与 ≥5 个其他蛋白的 UKB-PPP trans-pQTL 可信集共享变异的 cs（多效热点）；
7. PAV 替代定义（规范转录本；PIP≥0.2；先导 PAV）。

## 6. 应用示例（描述性，非假设检验）
对常用心血管代谢药物靶点（PCSK9、HMGCR、NPC1L1、APOB、LPA、ANGPTL3、ANGPTL4、CETP、APOC3、LPL、LDLR、IL6R、IL6、IL1B、F2、F10、F11、ACE、AGT、REN、DPP4、GCG、LEP、FGF21、GDF15、PLA2G7、MPO、CRP 等，以 UKB-PPP 有 cis cs 者为准）的全部 cis 工具变量给出分级，写成“如何挑选更可靠的工具变量”，不评价他人既往结论。

## 7. 多重比较、软件与可重复性
- Holm：确认族 6 个主终点；各假设次要终点族内 Holm；敏感性与应用不校正（描述性）。
- 事后或探索性分析一律标注“事后/探索性”。
- 软件：Python（duckdb、pandas、numpy、scipy、statsmodels、scikit-learn），版本在运行时写入结果工作簿；随机种子 20261103。
- 偏离：任何数据处理错误的更正以“修订记录”追加，注明是否影响结局；终点与阈值不变。

## 8. 未获得数据（写入论文局限）
deCODE 全量汇总统计（用户决定不申请）；ARIC/Fenland 补充表；Seer 全基因组汇总统计；eQTL Catalogue lbf_variable；GTEx 全名义统计；OT GWAS×QTL 共定位。共定位使用 OT 99% 可信集内逐变异 logBF（非全区域 lbf），属近似。

## 修订记录
（空）
