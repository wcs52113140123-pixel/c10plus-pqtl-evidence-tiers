# Changelog

## v1.3 (2026-11-04) — Zenodo doi:10.5281/zenodo.23943481
- **Resource v2.3** (`resource/C10plus_cis_pqtl_robustness_grades_v2.3.tsv.gz`, SHA-256 `067b7f7c44767b1befa429e5721dafac8e7d5ecd8a88473ffd2e7ce550b3c8be`): the 76 v2.2 columns are unchanged value for value; 25 columns from addendum A4 were appended. These are two mass-spectrometry layers that are independent of the Seer study: Generation Scotland (serum SWATH, N = 15,035) and HOLBAEK (plasma Orbitrap Astral, N ≈ 1,909).
- **Analysis plan addendum A4** (`prereg/C10plus_prereg_addendum_A4_independent_MS.md`, SHA-256 `cc2da0441cba82e1248ff6c328ac692cc5043fa0b6c133a5eb30f74c01935fd7`, fixed 2026-11-04 06:36:03 UTC) and its code (`src/analysis/A4/`).
- **Post hoc analyses L1/L2** (`src/posthoc_L1L2/`; not part of any fixed analysis plan).
- **Field dictionaries**: English translation added (`resource/FIELD_DICTIONARY_en.md`); Chinese dictionary updated to v2.3.
- **Lookup portal** (`docs/`): resource v2.3; cards show MS support for four layers (MSPA, Seer, Generation Scotland, HOLBAEK); rsID search from the Open Targets 26.09 variant index; About page updated.
- **Metadata**: new archive title (“resource, code and analysis plans”); README, CITATION.cff and DATA_LICENSE.md updated (Generation Scotland, HOLBAEK and Niu et al. 2025 attributions). Contact e-mail and ORCID entries were removed from repository files.
- Code paths in the added scripts are relative to the analysis working directory.

## v1.2 (2026-11-04) — Zenodo doi:10.5281/zenodo.23938291
- Resource v2.2 with addendum A2 (ARIC and Fenland, SomaScan v4); analysis and figure code; lookup portal.
