args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: Rscript infercnv_run.R counts.tsv annotations.tsv gene_order.tsv normal_groups output_dir")
library(infercnv)
groups <- strsplit(args[4], ",", fixed = TRUE)[[1]]
if (!length(groups) || any(!nzchar(groups))) stop("At least one normal group is required")
obj <- CreateInfercnvObject(raw_counts_matrix = args[1], annotations_file = args[2],
  delim = "\t", gene_order_file = args[3], ref_group_names = groups)
infercnv::run(obj, cutoff = 0.1, out_dir = args[5], cluster_by_groups = TRUE,
  denoise = TRUE, HMM = TRUE, HMM_type = "i6", num_threads = 4)
