# Figure code, release v1.3.3

This code regenerates the redesigned main figures (Fig. 1–6) and Supplementary Figs S1–S8 in `figures/`. It supersedes `../C10plus_figures_code_v1.2.py`.

```bash
pip install -r ../../../requirements.txt   # also needs scikit-learn and openpyxl
cp figure_inputs.example.json figure_inputs.json   # then edit the three paths
python run_all.py --config figure_inputs.json
# or: python run_all.py --workbook <xlsx> --resource <tsv.gz> --af2 <xlsx>
```

| Input | Source |
|---|---|
| `resource` | `resource/C10plus_cis_pqtl_robustness_grades_v2.3.tsv.gz` (this repository) |
| `workbook` | results workbook (Additional file 5 of the paper, v6) |
| `af2` | Additional file 2 of the paper; sheets Table S23–S36 are read as figure inputs |

Relative paths in the config file are resolved against the config file's folder. Output (PDF, PNG and TIFF per figure, plus `figure_layout_report.json`) goes to `out/` next to the code.

Modules: `c10style.py` (style and layout checks), `common.py` (shared helpers), `fig1.py`–`fig6.py`, `sfigs_a.py` (S1–S4), `sfigs_b.py` (S5–S8), `run_all.py` (driver).
