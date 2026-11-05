# Figures

`Fig4_mass_spec_four_layers.png` / `.pdf` (release v1.3.2) replaces the earlier `Fig4_mass_spec_truth_layer.png`. The figure covers all four mass-spectrometry (MS) layers: MSPA and the Seer full summary statistics (both from the Seer Proteograph study), plus Generation Scotland (serum SWATH) and HOLBAEK (plasma Orbitrap Astral), the two independent MS cohorts added in addendum A4 (resource v2.3).

- Platform concordance was positively associated with MS support in the independent cohorts: Generation Scotland OR 1.86, P = 0.009; HOLBAEK OR 4.67, P = 1.4 × 10⁻⁵.
- Transcript (eQTL) support did not predict MS confirmation.
- Signals whose 95% credible set contains a protein-altering variant (PAV) were confirmed less often in all four layers (independent cohorts combined: OR 0.61, P = 0.056).

The other PNGs (Fig. 1–3, 5, 6) were produced by `src/figures/C10plus_figures_code_v1.2.py` for release v1.2. Their MS panels, where present, use only the two Seer-study label sets.

In `src/figures/C10plus_figures_code_v1.2.py`, the two v1.2 Fig. 4 functions now write a superseded Seer-only panel (`FigS_v1.2_Seer_MS_layers_superseded`) with a neutral title. They no longer write a file named Fig. 4 or a title stating a conclusion.
