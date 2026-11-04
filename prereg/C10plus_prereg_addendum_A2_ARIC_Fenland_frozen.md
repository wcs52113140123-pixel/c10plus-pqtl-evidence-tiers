# C10+ 预注册增补 A2（ARIC + Fenland：独立 SomaScan v4 队列复现 H2）——冻结版
**冻结时间**：2026-11-04 04:32:55 UTC
**母预注册**：C10plus_preregistration_v2_frozen.md（SHA256 126f515f…62a89b0）；增补 A1（SHA256 6d54bd54…035e2）。除本增补明确定义者外，所有定义、协变量、模型、阈值沿用母预注册 H2。随机种子 20261103。

## English summary
Addendum A2 tests whether the confirmatory H2 results (deCODE SomaScan v4: eQTL anchoring predicts cross-affinity replication, H2a; protein-altering variants raise binary replication after power matching, H2b) reproduce in two further SomaScan v4 cohorts: ARIC European ancestry (Zhang 2022 Nat Genet, n≈7,213; sentinel, conditionally independent and fine-mapped cis-pQTLs, hg38) and Fenland (Pietzner 2021 Science, n=10,708; GWAS Catalog curated significant lead associations, GRCh38, protein names mapped to genes by a fixed rule). Replication of a UKB-PPP cis credible set = any cognate-gene cis signal variant from the cohort falls in the UKB 95% credible set (position match). Endpoints, wrong-gene negative control, power matching and thresholds are those of H2; the family {H2a-ARIC, H2b-ARIC, H2a-Fenland, H2b-Fenland} is Holm-adjusted at α=0.05. Secondary: ARIC African ancestry, effect-direction agreement (H2c analogue), and DerSimonian–Laird random-effects meta-analysis across deCODE, ARIC-EA and Fenland. **The deCODE H2 results are known (non-blind for the hypotheses); ARIC and Fenland outcomes were not computed before freezing.** All three cohorts use SomaScan v4, so A2 adds reproducibility across cohorts, not orthogonality across platforms.

## 0. 盲法与已知信息
1. **已知**：母预注册 H2（deCODE）结果——H2a OR 3.35（2.63–4.28），H2b 匹配后 OR 2.74（2.31–3.25），H2c PAV 效应相关更低、反向更多；以及 A1 全部结果。因此 A2 的方向预测来自 H2，对假设非盲。
2. **写定前只看了**：ARIC 补充表 MOESM4 的表名、表头与行数；Fenland 研究/关联表的列名与行数；以及按 §2 规则得到的映射计数与映射率（ARIC-EA：2,004 个有 cis 信号的 SOMAmer，2,003 个映射到基因，1,933 个基因，其中 895 个为 UKB-PPP cis 基因；ARIC-AA：1,618 个 SOMAmer，1,567 个基因，755 个为 UKB cis 基因；Fenland：3,892 个研究中 3,618 个映射（93.0%），2,870 条 cis 关联，1,491 个基因，其中 776 个为 UKB cis 基因）。**未计算任何 ARIC/Fenland 信号与 UKB 可信集的重合，也未把它们与锚定或 PAV 做交叉表。**
3. **平台**：deCODE、ARIC、Fenland 都用 SomaScan v4（适配体基本相同）。因此 A2 检验的是“同一亲和平台在独立队列中的可复现性”，不增加平台正交性；PAV 型表位效应若存在于 SomaScan 适配体，会在三个队列中同样出现。
4. **样本重叠**：ARIC（美国社区队列）和 Fenland（英国剑桥郡，独立招募）都不是 UK Biobank 参与者来源，已知无重叠；但无法在个体层面核实。ARIC 与 Fenland 之间也无重叠。
5. **数据局限**：ARIC 只公布有显著 cis 信号的 SOMAmer（ST3），没有 LD 代理变异；Fenland 只有 GWAS Catalog 收录的显著先导关联，没有全量汇总统计，也没有 SOMAmer/UniProt ID，需要按名称映射。

## 1. 分析集
- 单位：UKB-PPP EUR cis 可信集（同 H2）。
- 每个队列的分析集：同源基因在该队列中至少有 1 个 cis 信号的 UKB cs（即 has_cis=1；H2 已证明 has_cis=0 时不可能复现、logistic 中完全分离，限制在 has_cis=1 与原模型 MLE 等价）。
- 主分析使用全部符合条件的基因；敏感性：只用 H2 的 888 个基因（v4 有、Sun 2018 无）。

## 2. 映射与信号定义
- **ARIC**：SOMAmer→基因：UniProt（OT gene.uniprot_swissprot 唯一命中）→ 否则 Entrez 基因符号唯一命中。信号变异 = ST3.x 先导 + ST6.x 条件独立变异 + ST10.x 精细定位变异（全部为 cis，hg38）。
- **Fenland**：研究性状名标准化（小写、统一引号、去掉末尾 “levels”），依次匹配：(1) ARIC 表中 SomaScan 靶蛋白短名/全名（精确，其次去掉末尾括号）；(2) UniProt Swiss-Prot 人类推荐名与别名（精确，其次去掉末尾括号；仅取映射到唯一 Ensembl 基因的条目）。一名多基因判为歧义，不映射；复合物（含 “:”、heterodimer/trimer/tetramer、complex）除非经 (1) 映射，否则不映射。**不使用 GWAS Catalog 的 mapped_genes**（那是变异邻近基因，用于映射会循环）。cis = 先导位于同源基因 [start−1 Mb, end+1 Mb]（GRCh38，OT 坐标）。
- 映射率：见 §0.2（93.0%）。

## 3. 复现判定
- **ARIC 主规则** `rep_aric`：同基因任一 ARIC 信号变异（先导/条件独立/精细定位）的 hg38 位置落在 UKB 95% 集内。次要：(a) 严格——ST10 精细定位变异且 PIP≥0.5 落在 UKB 95% 集；(b) 放宽——任一信号变异落在 UKB 99% 集。
- **Fenland 主规则** `rep_fenland`：同基因任一 cis 先导变异落在 UKB 95% 集（GRCh38 位置；位置缺失时用 rsID 经 OT variant_core 映射）。次要：落在 99% 集。
- 注：两层都没有 LD 代理变异，只能做“精确变异重合”，比 deCODE 的 LD 类规则更严格；预期复现率更低，这是非差异性的保守偏差。

## 4. 终点（沿用 H2）
- **H2a-ARIC / H2a-Fenland**：logistic `rep ~ A_gtex + mlp + logn95 + PAV`（分析集内），基因聚类 SE；成立需 OR ≥ 1.5、Holm P < 0.05，且联合模型中错基因对照 `W_gtex` 的 OR 显著更小（单侧 P < 0.05）。次要：A_nongtex；功效匹配子集（UKB P<1e-30）OR ≥ 1.3。
- **H2b-ARIC / H2b-Fenland**：PAV 功效匹配——按 (mlp, absb, maf, logn95) 做 1:3 Mahalanobis 最近邻、无放回、标准化后距离 ≤ 0.5（与 H2 实现相同）；匹配样本上 logistic `rep ~ PAV`，按匹配集聚类；成立需 OR ≥ 1.3 且 Holm P < 0.05。同时报告样条调整 OR 与匹配前后 SMD。
- **增补族 Holm**：{H2a-ARIC, H2b-ARIC, H2a-Fenland, H2b-Fenland}，α=0.05。判定规则同母预注册。
- **停止规则**：某队列分析集中复现事件 < 30 → 该队列 H2a/H2b 判为不确定；H2b 匹配集 < 50 → H2b 不确定。

## 5. 次要（预设，不计入确认族）
1. **ARIC-AA**（跨祖先 + 跨平台）：用 ST3.2/6.2/10.2 重复 H2a、H2b。
2. **H2c 类效应一致性**：ARIC 用匹配变异的单 SNP β（测试等位基因）；Fenland 用先导 β 与效应等位基因；UKB 用该变异的 locus β（先导用 qtl_cs.beta）；等位基因协调：测试/效应等位基因 == alt → +，== ref → −，否则剔除。按 PAV 分层的 Pearson r（Fisher z，单侧预测 PAV 更低）与方向相反比例（Fisher，单侧预测 PAV 更高）。
3. **三队列随机效应 meta**（DerSimonian–Laird）：H2a 与 H2b 的 log OR，合并 deCODE（H2，已知）、ARIC-EA、Fenland；报告合并 OR、95% CI、I²；另报告不含 deCODE 的 ARIC-EA + Fenland 合并。
4. 敏感性：只用 H2 的 888 个基因；ARIC 严格/放宽规则；Fenland 99% 规则。

## 6. 报告方式
A2 作为预注册增补，与 H2 并列（同一表中 deCODE / ARIC-EA / Fenland / meta 各列）。论文中写明三者同为 SomaScan v4，A2 只证明“这一模式在同平台独立队列中可复现”，不能当作正交验证。

## 7. 资源表新增字段（v2.2）
aric_ea_has_cis、rep_aric_ea、rep_aric_ea_strict、dir_aric_ea、aric_aa_has_cis、rep_aric_aa、fenland_has_cis、rep_fenland、dir_fenland。

## 修订记录
（空）
