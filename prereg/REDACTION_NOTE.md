# Redaction note for the published analysis plans (release v1.3.2)

The analysis plans in this folder are **redacted copies** of the plans as they were fixed before the analyses they govern.

**What changed.** Only lines that named computing infrastructure were changed: machine names and aliases, remote-login aliases, a Linux user name and absolute paths on the analysis machine. They now read "a local Linux workstation", "a rented cloud server" or `<project_dir>/`. Three lines in three plans were changed (one line each). Plans A2 and A3 contained no such lines and are byte-identical to the frozen originals.

**What did not change.** Scientific content is unchanged: hypotheses, endpoints, definitions, covariates, models, thresholds, multiplicity rules, decision rules and all dates and freeze timestamps. Each changed line keeps its line number, so the line count of every file is unchanged.

**Checksums.** The SHA-256 values recorded at freeze (column "Original SHA-256") are the ones cited in the paper, in the other plans and in the resource README. They identify the unredacted originals. The per-plan `.sha256` files in this folder are those freeze records, unchanged. They therefore match the published files only for A2 and A3. To verify the files as published, use `SHA256SUMS_redacted.txt` (`sha256sum -c SHA256SUMS_redacted.txt`).

**Unredacted originals** are available to journal editors and reviewers on request from the corresponding author. Their SHA-256 values match the "Original SHA-256" column.

| Plan | Scope | Frozen at | Original SHA-256 (as recorded at freeze) | SHA-256 of redacted copy | Lines changed |
|---|---|---|---|---|---|
| `C10plus_preregistration_v2_frozen.md` | Main analysis plan v2.0 (H1–H5) | 2026-11-03 11:02:32 UTC | `126f515f4734baef319629584521bc492e4a70730067f293c76964c2b62a89b0` | `f4e33330285ee649c9a8c2c392601dde87a1c9e0a2fd316d49c66c4f0046995e` | 1 |
| `C10plus_prereg_addendum_A1_seer_full_frozen.md` | Addendum A1 (Seer full summary statistics) | 2026-11-03 12:25:37 UTC | `6d54bd5478c43184e349523b95c0406ded7bab5e106823e0d7537ad85b2035e2` | `ef85628d064be5f4c9c73a625132c6fe2ca89bd33386651ad9b084f7880ad423` | 1 |
| `C10plus_prereg_addendum_A2_ARIC_Fenland_frozen.md` | Addendum A2 (ARIC and Fenland, SomaScan v4) | 2026-11-04 04:32:55 UTC | `08df92f09ddaa5e82e54e87e134cab325dc6d02224b00456dcde6599522cbdb3` | `08df92f09ddaa5e82e54e87e134cab325dc6d02224b00456dcde6599522cbdb3` | 0 |
| `C10plus_prereg_addendum_A3_practical_impact.md` | Addendum A3 (practical impact on published MR instruments) | 2026-11-04 06:07:53 UTC | `4a5170ac229e012b253a36046f7163144e0669ec1fec8f66b4e77d81279c3971` | `4a5170ac229e012b253a36046f7163144e0669ec1fec8f66b4e77d81279c3971` | 0 |
| `C10plus_prereg_addendum_A4_independent_MS.md` | Addendum A4 (Generation Scotland and HOLBAEK MS layers) | 2026-11-04 06:36:03 UTC | `cc2da0441cba82e1248ff6c328ac692cc5043fa0b6c133a5eb30f74c01935fd7` | `c3578de084c70dc68b679e71352e34467ea3041637f2df57c926cb18f4358df7` | 1 |

A machine-readable, line-level summary of the changes is in `redaction_line_diff.json`. It gives the plan, line number and category of each change, never the removed text.
