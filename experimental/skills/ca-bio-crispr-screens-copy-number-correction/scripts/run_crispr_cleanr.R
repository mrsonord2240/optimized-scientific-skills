# Reference: CRISPRcleanR 3.0+, BiocManager | Verify API if version differs
# Function signatures checked against francescojm/CRISPRcleanR upstream source, 2026-09-16.
#
# CRISPRcleanR pre-hoc copy-number bias correction for a single cancer-line screen.
# Outputs corrected counts compatible with MAGeCK / BAGEL2 / drugZ downstream.

library(CRISPRcleanR)

# === INPUTS ===
counts_file <- 'counts.txt'                       # tab-separated: sgRNA, gene, sample columns
plasmid_col <- 'Plasmid'
sample_cols <- c('Veh_r1', 'Veh_r2', 'Drug_r1', 'Drug_r2')
data(KY_Library_v1.0)                              # built-in KY library; replace with custom

# === LOAD AND NORMALIZE ===
# ccr.NormfoldChanges takes a FILE PATH as its first argument (or NULL plus Dframe=)
norm <- ccr.NormfoldChanges(counts_file,
                              min_reads = 30,                # min reads/sgRNA: drop low-coverage
                              EXPname   = 'cancer_screen',
                              libraryAnnotation = KY_Library_v1.0)

# === GENOME-SORTED LFC ===
gw_lfc <- ccr.logFCs2chromPos(norm$logFCs,
                                KY_Library_v1.0)

# === CN-AWARE SEGMENTATION + CORRECTION ===
# ccr.GWclean detects spatial enrichment / depletion patterns and shifts segments
# toward the global mean. Unsupervised: no CN profile required.
cleaned <- ccr.GWclean(gw_lfc,
                         display = TRUE,
                         label   = 'cancer_screen')

# === DERIVE CORRECTED COUNTS FOR DOWNSTREAM TOOLS ===
# Positional signature: (CL, normalised_counts, correctedFCs_and_segments, libraryAnnotation, ...)
corrected_counts <- ccr.correctCounts('cancer_screen',
                                        norm$norm_counts,
                                        cleaned,
                                        KY_Library_v1.0,
                                        OutDir = './')
# Returns the corrected count data.frame and writes <CL>_correctedCounts.RData;
# write your own TSV for MAGeCK input.

# === DIAGNOSTIC: PRE vs POST CORRECTION ===
# ccr.logFCs2chromPos() returns CHR, startp, endp, genes, avgFC, BP -- there is no CN column and
# no column named gw_lfc, so a CN profile has to be merged in by gene.
# cn_profile.txt: two columns, gene and copy_number (from WGS / SNP-array / a cell-line database).

if (file.exists('cn_profile.txt')) {
  cn <- read.table('cn_profile.txt', header = TRUE, sep = '\t')   # columns: gene, copy_number

  pre  <- data.frame(gene = gw_lfc$genes,                 avgFC = gw_lfc$avgFC)
  post <- data.frame(gene = cleaned$corrected_logFCs$genes,
                     avgFC = cleaned$corrected_logFCs$correctedFC)

  pre  <- merge(pre,  cn, by = 'gene')
  post <- merge(post, cn, by = 'gene')
  amp_pre  <- pre$copy_number  > 4
  amp_post <- post$copy_number > 4
  stopifnot(sum(amp_pre) > 0)   # no amplified gene matched: check gene-name conventions

  cat('Pre-correction logFC mean (amplified):',  mean(pre$avgFC[amp_pre]),   '\n')
  cat('Post-correction logFC mean (amplified):', mean(post$avgFC[amp_post]), '\n')
  cat('Correction effectiveness: difference =',
      mean(post$avgFC[amp_post]) - mean(pre$avgFC[amp_pre]), '\n')
  cat('Spearman rho with CN, pre  =', cor(pre$avgFC,  pre$copy_number,  method = 'spearman'), '\n')
  cat('Spearman rho with CN, post =', cor(post$avgFC, post$copy_number, method = 'spearman'), '\n')
}

# === HIT CALLING ON CORRECTED COUNTS ===
# Feed cancer_screen_cleanr_corrected_counts.txt to MAGeCK / BAGEL2 / drugZ
# system('mageck test -k cancer_screen_cleanr_corrected_counts.txt ...')
