# Purpose: per-locus, per-gene effector-gene concordance scoring across the six core evidence
#          streams (fine-mapping, coloc, distance, PoPS, L2G, ABC/ENCODE-rE2G) at the canonical
#          thresholds, with a confidence tier per gene.
# Input:   TSV, one row per (locus, gene), columns: locus, gene, pip_top_variant,
#          credible_set_purity, coloc_pph4, distance_to_tss, pops_decile_rank, l2g_score,
#          abc_score, encode_re2g_score. Leave a cell NA when that stream is unavailable.
# Output:  TSV of candidates with concordance >= 3 (default: <input>.concordance.tsv), plus the
#          full scored table (<input>.scored.tsv) with an n_streams_available column.
# Usage:   Rscript concordance_scoring.R locus_candidates.tsv [out_prefix]
# Needs:   dplyr
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop('usage: Rscript concordance_scoring.R locus_candidates.tsv [out_prefix]')
in_path <- args[1]
out_prefix <- if (length(args) >= 2) args[2] else in_path

# Per-locus candidate gene table; one row per (locus, gene)
candidates <- read.table(in_path, header = TRUE, sep = '\t')

# Score each evidence stream against canonical thresholds. A stream with no data (NA) counts
# as not passed, so a missing stream lowers the tier instead of turning the whole row NA
# (which filter() would silently drop); n_streams_available reports which streams had data.
candidates <- candidates %>%
  mutate(
    pass_finemap = coalesce(pip_top_variant > 0.5 & credible_set_purity > 0.5, FALSE),
    pass_coloc = coalesce(coloc_pph4 >= 0.7, FALSE),
    pass_distance = coalesce(distance_to_tss <= 100000, FALSE),
    pass_pops = coalesce(pops_decile_rank == 1, FALSE),
    pass_l2g = coalesce(l2g_score >= 0.5, FALSE),
    pass_abc = coalesce(abc_score >= 0.02 | encode_re2g_score >= 0.5, FALSE),
    n_streams_available = as.integer(!is.na(pip_top_variant) & !is.na(credible_set_purity)) +
                          as.integer(!is.na(coloc_pph4)) + as.integer(!is.na(distance_to_tss)) +
                          as.integer(!is.na(pops_decile_rank)) + as.integer(!is.na(l2g_score)) +
                          as.integer(!is.na(abc_score) | !is.na(encode_re2g_score)),
    concordance = pass_finemap + pass_coloc + pass_distance +
                  pass_pops + pass_l2g + pass_abc,
    confidence_tier = case_when(
      concordance >= 5 ~ 'near_certain',
      concordance >= 4 ~ 'strong',
      concordance >= 3 ~ 'high',
      concordance >= 2 ~ 'suggestive',
      TRUE ~ 'associational_only'))

# Report
report <- candidates %>%
  filter(concordance >= 3) %>%
  arrange(desc(concordance), desc(l2g_score)) %>%
  select(locus, gene, concordance, n_streams_available, confidence_tier, l2g_score,
         pops_decile_rank, coloc_pph4)

write.table(candidates, paste0(out_prefix, '.scored.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(report, paste0(out_prefix, '.concordance.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
print(as.data.frame(report))
