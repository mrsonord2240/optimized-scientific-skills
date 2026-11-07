# flag_excluded_region.R -- programmatic MHC / chr 8 inversion gate for coloc.
#
# Purpose: prose alone is not a safeguard. Call this on the locus BEFORE coloc.abf or
#          coloc.susie; it returns 'MHC', 'chr8_inversion' or NA, and
#          stop_if_excluded_region() turns a non-NA flag into a hard stop().
# Inputs:  chr (6, '6' or 'chr6'), pos_bp (hg38 base-pair position), build (hg38 only).
# Usage:   source('scripts/flag_excluded_region.R')
#          stop_if_excluded_region(chr = gwas_df$CHR[1], pos_bp = gwas_df$POS[which.min(gwas_df$P)])
#          (coloc_abf.R and coloc_susie.R source this file and call the gate themselves.)

flag_excluded_region <- function(chr, pos_bp, build = 'hg38') {
  if (build != 'hg38') stop('flag_excluded_region: liftover to hg38 first -- this Skill only documents hg38 MHC / chr8 inversion boundaries')
  chr <- gsub('^chr', '', as.character(chr))
  if (chr == '6' && pos_bp >= 25000000 && pos_bp <= 35000000) return('MHC')
  if (chr == '8' && pos_bp >= 8100000  && pos_bp <= 11900000) return('chr8_inversion')
  NA_character_
}

stop_if_excluded_region <- function(chr, pos_bp, build = 'hg38') {
  region_flag <- flag_excluded_region(chr = chr, pos_bp = pos_bp, build = build)
  if (!is.na(region_flag)) {
    stop(sprintf(
      'Locus is in the %s exclusion zone -- standard coloc PP.H4 is not interpretable here. Use HLA-coloc (Butler-Laporte 2024) for MHC, or pre-condition on inversion genotype for chr8; do not report a naive coloc.abf/coloc.susie result.',
      region_flag))
  }
  invisible(NA_character_)
}
