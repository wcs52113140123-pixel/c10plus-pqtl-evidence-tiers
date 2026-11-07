# Figures (release v1.3.4)

These are the redesigned main and supplementary figures. Each is provided as vector PDF and 300-dpi PNG. They were produced by the code in `src/figures/v1.3.3/`; see its README for inputs and usage.

## Main figures

| File | Content |
|---|---|
| `Fig1.pdf` / `.png` | Concept schematic (a: abundance effect vs PAV binding artefact; not data), four-step study workflow (b) and all pre-specified endpoints with thresholds and verdicts (c). Drawn by `src/figures/v1.3.3/fig1v3.py`. |
| `Fig2.pdf` / `.png` | SomaScan replication by PAV and eQTL-anchor group in deCODE and Sun 2018; power matching; anchoring and wrong-gene control. |
| `Fig3.pdf` / `.png` | Cross-platform effect-size agreement and sign reversal for PAV+ and PAV− pairs, reproduced across the SomaScan v4 cohorts. |
| `Fig4.pdf` / `.png` | MS support across all four MS layers (MSPA, Seer, Generation Scotland, HOLBAEK): label agreement, concordance, PAV, anchoring, and unanchored-PAV contrasts. |
| `Fig5.pdf` / `.png` | Replication in non-European ancestries, grading-model AUC, and MS support or replication by tier. |
| `Fig6.pdf` / `.png` | Use of Tier 4 instruments in published drug-target MR studies, enrichment among significant results, and Tier 1–2 alternatives. |

## Supplementary figures (`supplementary/`)

| File | Content |
|---|---|
| `FigS1` | Pre-specified sensitivity analyses for five endpoints. |
| `FigS2` | Extraction of the Seer full summary statistics; MS label outcomes and MSPA–Seer agreement. |
| `FigS3` | MS support across the full concordance × anchor × PAV grid, per MS layer. |
| `FigS4` | Cross-ancestry replication by ancestry; effect-size slopes. |
| `FigS5` | Grading-model coefficients and transfer across targets. |
| `FigS6` | Resource overview: tier distribution, tier by SomaScan concordance, signal strength by tier. |
| `FigS7` | Anchoring sensitivity, plus PAV and concordance contrasts, in the Seer-study MS layers. |
| `FigS8` | PAV-aware tiers for 42 cardiometabolic drug targets. |

## Superseded figure code

`src/figures/C10plus_figures_code_v1.2.py` (with `figcommon.py`) produced the v1.2 figures, which were removed in v1.3.3. It is **superseded** by `src/figures/v1.3.3/` and is kept only for the record.
