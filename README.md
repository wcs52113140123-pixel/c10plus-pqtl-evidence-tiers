# PAV-aware evidence tiers for UKB-PPP cis-pQTL credible sets (v1.3.4)

**Online lookup portal:** https://wcs52113140123-pixel.github.io/c10plus-pqtl-evidence-tiers/

**Paper (in submission):** *Protein-altering variants inflate cross-platform replication of plasma protein QTLs*. Jiayao Ni, Wanqiang Zhang, Tao Zhang.
**Archive:** v1.3.4 https://doi.org/10.5281/zenodo.24041516 · all versions https://doi.org/10.5281/zenodo.23938290 · **Code:** https://github.com/wcs52113140123-pixel/c10plus-pqtl-evidence-tiers
**Licences:** code MIT (`LICENSE`); data and resource CC BY 4.0 (`DATA_LICENSE.md`).
**How to cite:** Ni J, Zhang W, Zhang T. PAV-aware evidence tiers for UK Biobank Pharma Proteomics Project plasma cis-pQTLs: resource, code and analysis plans (v1.3.4). Zenodo; 2026. https://doi.org/10.5281/zenodo.24041516. Please also cite the original data sources (`DATA_LICENSE.md`).

**Key message.** Protein-altering variants (PAVs) inflate binary Olink–SomaScan concordance, reproducibly in three SomaScan v4 cohorts (deCODE, ARIC, Fenland; one platform). Concordance is therefore not orthogonal validation for PAV-containing signals. Overall, platform concordance was positively associated with MS support in the two independent MS cohorts (Generation Scotland, HOLBAEK; combined OR 2.92) but only weakly, and not significantly, in the two Seer layers. Tiers annotate PAV and transcript (eQTL) support; they are not a validated reliability ranking, because transcript support did not predict mass-spectrometry confirmation (see `resource/README.md`).

**Data sources and terms.** Only derived per-credible-set statistics are released; no third-party summary-statistic files are redistributed. Sources: UKB-PPP (Open Targets 26.09, CC0; non-European data via Synapse syn51364943, registration required), Sun et al. 2018 and eQTL Catalogue (CC BY 4.0), GTEx v10, deCODE (Eldjarn et al. 2023), Suhre et al. 2025, Seer Proteograph (GWAS Catalog, CC0), ARIC (Zhang et al. 2022), Fenland (GWAS Catalog), Generation Scotland (Edinburgh DataShare, CC BY 4.0) and HOLBAEK (GWAS Catalog, CC0; Niu et al. 2025, CC BY 4.0). See `DATA_LICENSE.md` for the full table and required citations.

## Layout
- `resource/` — evidence-tier table v2.3 (TSV.GZ; Parquet in the Zenodo archive; v2.2 kept for reference), README with field definitions, field dictionaries (`FIELD_DICTIONARY_en.md`, `FIELD_DICTIONARY_zh.md`), `example_select_instruments.py|R`
- `prereg/` — analysis plans fixed before the corresponding analyses: plan v2.0 and addenda A1, A2, A3 and A4, each with its SHA-256 file as recorded at freeze. The published copies have computing-infrastructure lines redacted; scientific content is unchanged (see `prereg/REDACTION_NOTE.md`; `SHA256SUMS_redacted.txt` verifies the published files).
- `src/analysis/` — analysis pipeline (`RUN_ORDER.md`); `src/analysis/A4/` — addendum A4 (Generation Scotland and HOLBAEK MS layers; `README_A4.md`)
- `src/posthoc_L1L2/` — post hoc analyses (not part of any fixed analysis plan): L1, main and MS endpoints with and without chrX; L2, diagnostics of the H4 HOLBAEK result (chromosome-blocked cross-validation and leave-one-MS-layer-out transfer)
- `src/figures/v1.3.3/` — figure code for the v1.3.3 figure set (`run_all.py`; input paths via `--config` or CLI arguments); `src/figures/C10plus_figures_code_v1.2.py` is superseded
- `src/revision/` — post hoc robustness analyses (reproduces Additional file 6 from the resource)
- `figures/` — main figures Fig. 1–6 and `figures/supplementary/` Figs S1–S8 (PDF and PNG, release v1.3.4; see `figures/README.md`)
- `docs/` — lookup portal (GitHub Pages) and `DATA_SOURCES.md` (upstream data and retrieval notes)
- `CHANGELOG.md` — release history

Script paths are relative to the analysis working directory (`./`), with the raw database expected under `../data_v1/`; adjust them to your layout.

## Quick start
```bash
pip install -r requirements.txt
cd resource && python example_select_instruments.py PCSK9
```
