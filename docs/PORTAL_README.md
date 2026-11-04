# C10+ query portal (static)

`index.html`, `assets/` and `data/` form a fully static, dependency-free query page for the
C10+ PAV-aware cis-pQTL evidence-tier resource v2.2 (14,472 UKB-PPP Olink cis credible sets).

- No external scripts, fonts, CDNs, analytics, cookies or local storage (enforced by a
  Content-Security-Policy meta tag: `default-src 'self'`).
- `data/c10plus_tiers_v2.2.json` is generated from
  `resource/C10plus_cis_pqtl_robustness_grades_v2.2.tsv.gz` by `tools/build_portal_data.py`
  (display fields only, rounded; rounding rules and source SHA-256 are in `meta`).
- `.nojekyll` disables Jekyll so files are served as-is.
- Serve locally: `python -m http.server -d docs 8000` then open http://localhost:8000/.

Data: CC BY 4.0 (doi:10.5281/zenodo.23938291). Code: MIT.
