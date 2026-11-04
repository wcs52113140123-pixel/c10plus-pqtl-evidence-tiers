# PAV-aware evidence tiers for UKB-PPP cis-pQTL credible sets (v1.2)

**Online lookup portal:** https://wcs52113140123-pixel.github.io/c10plus-pqtl-evidence-tiers/

**Paper (in submission):** *Protein-altering variants inflate agreement between affinity proteomics platforms for plasma cis-pQTLs: a pre-registered evaluation against transcript, mass-spectrometry and cross-ancestry evidence*. Jiayao Ni, Wanqiang Zhang, Tao Zhang.
**Archive DOI:** https://doi.org/10.5281/zenodo.23938291 · **Code:** https://github.com/wcs52113140123-pixel/c10plus-pqtl-evidence-tiers
**Licences:** code MIT (`LICENSE`); data and resource CC BY 4.0 (`DATA_LICENSE.md`).
**How to cite:** Ni J, Zhang W, Zhang T. C10+ PAV-aware evidence tiers for UKB-PPP cis-pQTL credible sets (v1.2). Zenodo. https://doi.org/10.5281/zenodo.23938291. Please also cite the original data sources (`DATA_LICENSE.md`).
**Contact:** Tao Zhang, zht779100@njmu.edu.cn

**Key message.** Protein-altering variants inflate binary Olink–SomaScan concordance, reproducibly in three SomaScan v4 cohorts (deCODE, ARIC, Fenland; one platform). Concordance is therefore not orthogonal validation. Tiers are categorical evidence patterns, not a validity ranking (see `resource/README.md`).

**Data sources and terms.** Only derived per-credible-set statistics are released; no third-party summary statistics are redistributed. Sources: UKB-PPP (Open Targets 26.09, CC0; non-European data via Synapse syn51364943, registration required), Sun et al. 2018 and eQTL Catalogue (CC BY 4.0), GTEx v10, deCODE (Eldjarn et al. 2023), Suhre et al. 2025, Seer Proteograph (GWAS Catalog), ARIC (Zhang et al. 2022) and Fenland (GWAS Catalog). See `DATA_LICENSE.md` for the full table and required citations.

## Layout
- `resource/` — evidence-tier table v2.2 (TSV.GZ; Parquet in the Zenodo archive), README with field definitions, Chinese field dictionary, `example_select_instruments.py|R`
- `prereg/` — frozen pre-registration v2.0 and addenda A1 and A2, each with its SHA-256 file
- `src/analysis/` — analysis pipeline (`RUN_ORDER.md`)
- `src/figures/` — figure code (v1.2) and shared plotting helpers
- `src/revision/` — post hoc robustness analyses (reproduces Additional file 6 from the resource)
- `figures/` — main figures (PNG)
- `docs/DATA_SOURCES.md` — upstream data and retrieval notes

## Quick start
```bash
pip install -r requirements.txt
cd resource && python example_select_instruments.py PCSK9
```
