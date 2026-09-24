# Purpose: focused deterministic checks for raincloud_phd.R.
# Inputs: optional path to raincloud_phd.R. Usage: Rscript test_distribution_plots.R [example-path]

args <- commandArgs(trailingOnly = TRUE)
example_path <- if (length(args)) args[[1L]] else 'raincloud_phd.R'
example_path <- normalizePath(example_path, mustWork = TRUE)
out_dir <- tempfile('distribution-plots-test-')
dir.create(out_dir)
old_wd <- setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

source(example_path)

stopifnot(
  identical(safe_violin_bw(data.frame(group = c('A', 'A', 'B', 'B'), value = c(0, 0, 1, 2))), 'nrd0'),
  identical(levels(df_paired$cluster), c('Cluster 1', 'Cluster 2', 'Cluster 3')),
  all(table(df_paired$cluster, df_paired$condition) == 40L),
  identical(nonmissing_n_labels(data.frame(group = c('A', 'A', 'B'), value = c(1, NA, 2))),
            c('A\n(n=1)', 'B\n(n=1)')),
  all(unname(tapply(!is.na(df_med$value), df_med$group, sum)) == c(80L, 80L)),
  file.exists('raincloud.pdf'), file.info('raincloud.pdf')$size > 2000L
)

tied_panel <- data.frame(group = rep(c('Constant', 'Variable'), each = 12L),
                         value = c(rep(0, 12L), seq(-1, 1, length.out = 12L)))
tied_bw <- safe_violin_bw(tied_panel)
tied_plot <- ggplot(tied_panel, aes(group, value, fill = group)) +
  geom_violin(bw = tied_bw, trim = TRUE) + geom_point(position = position_jitter(seed = 1))
tied_build <- ggplot_build(tied_plot)
ggsave('tied-panel.png', tied_plot, width = 5, height = 3, dpi = 120)
stopifnot(identical(tied_bw, 'nrd0'), nrow(tied_build$data[[1L]]) > 0L,
          file.info('tied-panel.png')$size > 2000L)

annotation_plot <- ggplot(data.frame(group = c('A', 'A', 'B'), value = c(1, NA, 3)),
                          aes(group, value)) +
  stat_summary(geom = 'text', vjust = -0.4,
    fun.data = function(x) data.frame(y = max(x, na.rm = TRUE),
                                      label = paste0('n=', sum(!is.na(x)))))
annotation_data <- ggplot_build(annotation_plot)$data[[1L]]
stopifnot(identical(sort(annotation_data$label), c('n=1', 'n=1')),
          all(is.finite(annotation_data$y)))

rain_build <- ggplot_build(p_raincloud)
box_layer <- rain_build$data[[2L]]
expected_medians <- as.numeric(tapply(df_med$value, df_med$group, median, na.rm = TRUE))
stopifnot(isTRUE(all.equal(sort(box_layer$middle), sort(expected_medians), tolerance = 1e-10)))

ggsave('split_violin.png', p_split, width = 6, height = 4, dpi = 120)
stopifnot(file.info('split_violin.png')$size > 2000L)
cat('PASS: ggdist raincloud, N labels, safe bandwidth fallback, and guarded split violin\n')
