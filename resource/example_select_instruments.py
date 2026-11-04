"""Minimal example: choose cis-pQTL instruments by PAV-aware tier.
Requires pandas. Run: python example_select_instruments.py PCSK9
"""
import sys
import pandas as pd

gene = sys.argv[1] if len(sys.argv) > 1 else "PCSK9"
grades = pd.read_csv("C10plus_cis_pqtl_robustness_grades_v2.3.tsv.gz", sep="\t", low_memory=False)

# 1. All cis credible sets for the target, strongest first
cs = grades[grades.symbol == gene].assign(abs_beta=lambda d: d.lead_beta.abs())
cs = cs.sort_values(["lead_mlog10p", "abs_beta"], ascending=False)
n_tied = (cs.lead_mlog10p == cs.lead_mlog10p.max()).sum() if len(cs) else 0
if n_tied > 1:
    print(f"WARNING: {n_tied} credible sets tied at the top -log10 P; lead chosen by larger |beta| - inspect them individually")
print(cs[["studyLocusId", "lead_variant_id", "lead_mlog10p", "tier", "tier_label",
          "platform_concordance", "ms_primary", "ms_ext_primary", "xa_replicated",
          "rep_aric_ea", "rep_fenland", "flag_concordant_unanchored_PAV"]].to_string(index=False))

# 2. Strongest Tier 1 credible set (anchored, no PAV): candidate sensitivity instrument
#    (tiers are evidence patterns, not a validity ranking)
t1 = cs[cs.tier == 1]
print("\nStrongest Tier 1 credible set:",
      t1.iloc[0].lead_variant_id if len(t1) else "none - consider MS/functional evidence")

# 3. Genome-wide: concordant but unanchored PAV signals (key caution class)
flag = grades[grades.flag_concordant_unanchored_PAV == 1]
print(f"\n{len(flag)} credible sets in {flag.gene_id.nunique()} genes are flagged concordant_unanchored_PAV")

# 4. Exclude post hoc beta-inconsistent credible sets before effect-size work
clean = grades[grades.flag_beta_inconsistent == 0]
print(f"{len(clean)} of {len(grades)} credible sets pass the beta-consistency flag")
