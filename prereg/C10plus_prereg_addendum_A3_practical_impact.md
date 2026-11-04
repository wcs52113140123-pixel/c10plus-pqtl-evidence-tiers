# C10+ 预注册附录 A3：已发表药靶/蛋白组 MR 工具变量的实际影响分析

版本：A3 v1.0（冻结时间与 SHA256 见同名 .sha256 文件）。
上位文件：预注册 v2.0（SHA256 126f515f…62a89b0）；附录 A1、A2。
资源版本：C10plus_cis_pqtl_robustness_grades_v2.2（Tier 1–4 + flag_concordant_unanchored_PAV + flag_beta_inconsistent），冻结后不得修改资源分级规则。

## 0. 冻结时状态声明
冻结时已完成：文献检索、补充材料下载、表头检查、纳入论文名单的确定。
冻结时**未进行**：任何已发表工具变量与资源分级（Tier/flag）的匹配；未查看任何已发表工具所在 CS 的 Tier。

## 1. 文献纳入（已在冻结前确定，冻结后不再增删）
检索：Europe PMC REST，OPEN_ACCESS，2023–2026，三组检索式（UKB-PPP / Olink pQTL + Mendelian randomization / drug target），得 493 篇候选；标题含 CKM 结局（心血管/肾/代谢）+ MR 关键词筛至 91 篇，全部下载 Europe PMC 合法开放补充材料。
纳入标准：(a) 暴露为 UKB-PPP Olink cis（或全基因组）pQTL 工具；(b) 补充材料同时提供逐蛋白工具变量 rsID 与蛋白–结局 MR 结果（或可明确区分显著/不显著）；(c) 可区分 UKB-PPP 来源与其他平台来源。
纳入名单（12 个来源）：
| # | 来源 | DOI | PMID | 结局 | 工具表 | 结果表 | 论文自报显著性定义 |
|---|---|---|---|---|---|---|---|
| S0 | Sun 2023 Nature ST30（UKB-PPP 原文 MR 部分） | 10.1038/s41586-023-06592-6 | 37794186 | 血脂/MI/CHD/卒中 | ST30 rs11591147 | ST30 | MR P<0.05 |
| 1 | JVS Vasc Sci 2025 AAA | 10.1016/j.jvssci.2025.100290 | 40607238 | AAA | Table_S1 study=ukb | Table_S2 dset=UKB | padj<0.05 |
| 2 | Nat Commun 2025 AF | 10.1038/s41467-025-61720-2 | 40645996 | AF | Table S20 Source=UKBPPP | Table S10 SMR | P_Bonferroni<0.05 |
| 3 | Mol Genet Genomics 2026 跨人群 pwMR | 10.1007/s00438-026-02458-4 | 42545500 | AF/CAD/HF/IHD/PAD（仅 European） | MOESM1 *_European | MOESM3 Table S3-1（显著列表；未列出者视为不显著） | FDR<0.05 |
| 4 | BMC Med 2026 空气污染→CKD | 10.1186/s12916-026-04909-6 | 42106700 | CKD | Table S26 | Table S27 IVW/Wald | FDR<0.05 |
| 5 | Biol Sex Differ 2026 血脂 | 10.1186/s13293-026-00934-5 | 42304513 | HDL/LDL/TC/TG/nonHDL | ST1 IV_list.all | ST1 *.all | fdr.all<0.05 |
| 6 | J Transl Med 2024 CHD/MI | 10.1186/s12967-024-05178-8 | 38689297 | CHD、MI | Table S19 Author=Sun_2 | Table S1（CHD）/S3（MI）Author=Sun_2 | P<0.05/检验数（Bonferroni，按论文） |
| 7 | Cardiovasc Diabetol 2025 PWAS-SMR | 10.1186/s12933-025-02847-w | 40751156 | 多种 CVD | Supp Table 6 | Supp Table 7 | SMR_causal==YES（若不存在 YES 则 SMR.multi.P<0.05/行数） |
| 8 | Cardiovasc Ther 2026 冠脉粥样硬化 | 10.1155/cdr/7359278 | 42666106 | 冠脉粥样硬化 | Table 1 | Table 2–4（显著列表；未列出者视为不显著） | 论文主分析显著 |
| 9 | Genome Med 2024 CKD | 10.1186/s13073-024-01356-x | 38898508 | CKD/eGFR 相关 | Table S6 data=UK Biobank | Table S7 data=UK Biobank | P<0.05 |
| 10 | Diabetes Obes Metab 2026 T2D 聚类 | 10.1111/dom.71253 | 42638160 | CAD/MI/AF/HF 等 | ST21 | ST23 | P<0.05 |
| 11 | EBioMedicine 2025 心代谢蛋白–骨骼肌肉 | 10.1016/j.ebiom.2025.105580 | 39919333 | 骨/肌肉结局（非 CKM，列为 other） | cis SNPs 表（Sun et al.） | 对应 MR 表 | 论文报告的 FDR/P 阈值 |
排除（冻结前确定）：非 UKB-PPP 暴露（deCODE/ARIC/SCALLOP 等）或未提供逐蛋白工具 rsID 的论文（如 PMC12983497、PMC12529848、PMC12124044 仅 10 蛋白且含 trans 工具、PMC13099373 非 pQTL 工具）。

## 2. 单位与协调
- 分析单位：已发表的 (来源, 蛋白, 结局) 对（"protein–outcome pair"）；工具单位：(来源, 蛋白, rsID)。
- 多方法时取论文主方法（IVW；单工具取 Wald ratio；SMR 论文取 SMR）。多祖先仅取欧洲人群。性别分层论文取合并性别结果。
- 显著性一律采用论文自报阈值（上表），不重新计算。
- 蛋白→基因：优先 UniProt，其次 Olink OID，再次 HGNC 符号（经 OT 26.09 基因表，含同义名）。

## 3. 匹配规则（冻结）
- 变体定位：rsID → GRCh38 位置（OT 26.09 variant_core rsIds）；rsID 不可解析时用论文报告的 chr:pos（hg19 时用 UCSC hg19→hg38 chain liftover；hg38 直接用），并要求等位基因一致或互补。
- 主匹配：工具变体位于同一基因（任一 UKB-PPP Olink assay）的 UKB-PPP 欧洲 cis 95% CS 内；若落入多个 CS，取该变体 PIP 最高者所在 CS。
- 次匹配（敏感性）：落入 99% CS；或与 CS lead 同位置。无 LD 参考面板，LD 匹配不作主分析（偏离父指令中"位置/LD"之 LD 部分，预先声明）。
- 未匹配工具（trans、cis 但不在任何 CS 内）单独计数报告，不纳入分母。
- 对层面：对中≥1 个工具主匹配成功即为"已匹配对"。对暴露定义：Tier4_any = 任一匹配工具所在 CS 为 Tier 4；flag_any = 任一匹配 CS flag_concordant_unanchored_PAV=True；另报告 best_tier（最佳）与 worst_tier（最差）。

## 4. 主终点
- P1（描述性主终点）：已匹配对中 Tier4_any 比例与 flag_any 比例，Wilson 95% CI；同时报告工具层面比例。预设解读：Wilson 下限 ≥5% 判为"实质性比例"；<5% 判为"影响有限"。
- P2（检验性主终点，H-A3b）：显著 vs 不显著对中 Tier4_any 比例差异。主模型：logistic(Tier4_any ~ significant + C(来源))，按基因聚类稳健 SE；论文内无变异的来源自然不贡献信息。敏感性：合并 Fisher 精确检验。同样对 flag_any 检验。2 个检验 Holm 校正（α=0.05，双侧）。预设方向期望：显著对中更高（OR>1）；Holm P<0.05 且 OR>1 → 支持；Holm P<0.05 且 OR<1 → 方向相反；否则不支持。若 logistic 因分离/不收敛失败，以 Fisher 为主（预先声明）。

## 5. 次要终点（全部探索性，不做多重校正结论）
- 按结局类别（心血管 / 肾 / 代谢-血脂 / 其他）、发表年份（2023–24 / 2025 / 2026）、是否有 SomaScan 交叉验证（来源层面：1,2,6,7,9,10 有；其余无）分层的 P1 比例。
- 对 Tier4_any 或 flag_any 的对：同一蛋白是否存在 Tier 1–2 的替代 CS（alt_T12_available），比例与 CI。
- 单工具（Wald ratio）对子集重复 P1/P2。
- 与资源整体 Tier 4 基线比例（全部 UKB-PPP cis CS 中 Tier 4 比例）比较（二项检验）。
- 列出显著且 Tier4_any 的药靶对（应用实例表）。

## 6. 停止规则
已匹配对 <200 → A3 判为"不可下结论"，仅作描述性报告。

## 7. 解释边界
Tier 4 表示"仅凭亲和平台无法区分丰度效应与表位效应"，不等于工具无效或 MR 结论错误；结论措辞为"帮助选择更可靠的工具变量"，不对任何已发表研究作否定判断。
