# Minimal example: choose cis-pQTL instruments by PAV-aware tier.
# Requires readr and dplyr. Run: Rscript example_select_instruments.R PCSK9
suppressPackageStartupMessages({library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
gene <- if (length(args) > 0) args[1] else "PCSK9"
grades <- read_tsv("C10plus_cis_pqtl_robustness_grades_v2.2.tsv.gz", show_col_types = FALSE, guess_max = 20000)

cs <- grades %>% filter(symbol == gene) %>% arrange(desc(lead_mlog10p), desc(abs(lead_beta)))
n_tied <- sum(cs$lead_mlog10p == max(cs$lead_mlog10p))
if (n_tied > 1) cat("WARNING:", n_tied, "credible sets tied at the top -log10 P; lead chosen by larger |beta| - inspect them individually\n")
print(cs %>% select(studyLocusId, lead_variant_id, lead_mlog10p, tier,
                    platform_concordance, ms_primary, ms_ext_primary, xa_replicated,
                    rep_aric_ea, rep_fenland, flag_concordant_unanchored_PAV), width = Inf)

t1 <- cs %>% filter(tier == 1)
# Strongest Tier 1 credible set: candidate sensitivity instrument (tiers are evidence patterns, not a validity ranking)
cat("\nStrongest Tier 1 credible set:",
    if (nrow(t1) > 0) t1$lead_variant_id[1] else "none - consider MS/functional evidence", "\n")

# Tier distribution by SomaScan concordance
print(grades %>% count(platform_concordance, tier) %>%
        group_by(platform_concordance) %>% mutate(frac = round(n / sum(n), 3)))
