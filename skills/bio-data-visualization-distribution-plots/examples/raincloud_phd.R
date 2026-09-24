# Purpose: standalone, deterministic distribution-plot examples for small, medium,
# large, and paired two-condition groups. Inputs: none (replace the generated frames
# with project data that use the documented columns). Usage: Rscript raincloud_phd.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggbeeswarm)
  library(ggdist)
  library(lvplot)
  library(introdataviz)
})

okabe_subset <- c('#0072B2', '#D55E00', '#009E73')

nonmissing_n_labels <- function(data, group = 'group', value = 'value') {
  counts <- tapply(!is.na(data[[value]]), data[[group]], sum)
  sprintf('%s\n(n=%d)', names(counts), as.integer(counts))
}

safe_violin_bw <- function(data, group = 'group', value = 'value') {
  values <- split(data[[value]], data[[group]], drop = TRUE)
  sj_ok <- vapply(values, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L || length(unique(x)) < 2L) return(FALSE)
    bw <- tryCatch(stats::bw.SJ(x), error = function(e) NA_real_)
    is.finite(bw) && bw > 0
  }, logical(1))
  if (all(sj_ok)) 'SJ' else 'nrd0'
}

prepare_split_violin_data <- function(data, cluster = 'cluster',
                                      condition = 'condition', value = 'expression',
                                      min_n = 30L) {
  stopifnot(is.data.frame(data), all(c(cluster, condition, value) %in% names(data)))
  condition_values <- data[[condition]]
  if (!is.ordered(condition_values) || nlevels(condition_values) != 2L) {
    stop('condition must be an ordered factor with exactly two levels', call. = FALSE)
  }
  keep_rows <- !is.na(data[[cluster]]) & !is.na(condition_values) & is.finite(data[[value]])
  filtered <- data[keep_rows, , drop = FALSE]
  counts <- xtabs(rep(1L, nrow(filtered)) ~ filtered[[cluster]] + filtered[[condition]])
  keep_clusters <- rownames(counts)[apply(counts >= min_n, 1L, all)]
  filtered <- filtered[filtered[[cluster]] %in% keep_clusters, , drop = FALSE]
  filtered[[cluster]] <- droplevels(factor(filtered[[cluster]], levels = keep_clusters))
  filtered[[condition]] <- ordered(filtered[[condition]], levels = levels(condition_values))
  filtered
}

set.seed(20260923)
df_small <- data.frame(
  group = factor(rep(c('Ctrl', 'Low', 'High'), each = 15L), levels = c('Ctrl', 'Low', 'High')),
  value = c(rnorm(15, 4.5, 0.7), rnorm(15, 5.2, 0.5), rnorm(15, 6.1, 0.8))
)
df_med <- data.frame(
  group = factor(rep(c('Control', 'Treated'), each = 80L), levels = c('Control', 'Treated')),
  value = c(rnorm(80, 5.0, 0.7), rnorm(40, 3.2, 0.45), rnorm(40, 6.9, 0.55))
)
df_large <- data.frame(
  group = factor(rep(c('Ctrl', 'Low', 'High'), each = 250L), levels = c('Ctrl', 'Low', 'High')),
  value = c(rnorm(250, 3, 1), rnorm(250, 5, 1.1), rnorm(250, 4, 1.4))
)

split_cell <- function(cluster, condition, n, mean) {
  data.frame(cluster = cluster, condition = condition, expression = rnorm(n, mean, 0.8))
}
df_paired_raw <- rbind(
  split_cell('Cluster 1', 'Control', 40, 3.5), split_cell('Cluster 1', 'Treatment', 40, 4.6),
  split_cell('Cluster 2', 'Control', 40, 4.0), split_cell('Cluster 2', 'Treatment', 40, 5.4),
  split_cell('Cluster 3', 'Control', 40, 5.2), split_cell('Cluster 3', 'Treatment', 40, 5.9),
  split_cell('Missing treatment', 'Control', 40, 4.2),
  split_cell('Too small', 'Control', 20, 4.1), split_cell('Too small', 'Treatment', 20, 4.9)
)
df_paired_raw$cluster <- factor(df_paired_raw$cluster)
df_paired_raw$condition <- ordered(df_paired_raw$condition, levels = c('Control', 'Treatment'))
df_paired <- prepare_split_violin_data(df_paired_raw, min_n = 30L)

# 1. SMALL N (10-29): every point plus the group median.
small_labels <- nonmissing_n_labels(df_small)
p_small <- ggplot(df_small, aes(group, value, color = group)) +
  geom_quasirandom(method = 'quasirandom', width = 0.3, alpha = 0.7, size = 1.5) +
  stat_summary(fun = median, geom = 'crossbar', width = 0.5, color = 'black', linewidth = 0.4) +
  scale_color_manual(values = okabe_subset) +
  scale_x_discrete(labels = small_labels) +
  labs(x = NULL, y = 'Expression') +
  theme_classic(base_size = 10) + theme(legend.position = 'none')

# 2. MEDIUM N (30-200): ggdist works on current ggplot2; raw points are seeded.
medium_labels <- nonmissing_n_labels(df_med)
p_raincloud <- ggplot(df_med, aes(group, value, fill = group, color = group)) +
  stat_halfeye(adjust = 1, width = 0.6, .width = 0, justification = -0.2,
               point_colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = 'black') +
  geom_point(aes(x = as.numeric(group) - 0.2),
             position = position_jitter(width = 0.05, height = 0, seed = 20260923),
             alpha = 0.5, size = 1.2, shape = 16) +
  scale_fill_manual(values = c(Control = '#0072B2', Treated = '#D55E00')) +
  scale_color_manual(values = c(Control = '#0072B2', Treated = '#D55E00')) +
  scale_x_discrete(labels = medium_labels) +
  coord_flip() + labs(x = NULL, y = 'Biomarker level') +
  theme_classic(base_size = 10) + theme(legend.position = 'none')

# 3. LARGE N (>=201): letter-value plot retains tail structure.
large_labels <- nonmissing_n_labels(df_large)
p_lv <- ggplot(df_large, aes(group, value, fill = group)) +
  geom_lv(k = 5, alpha = 0.7) +
  scale_fill_manual(values = okabe_subset) +
  scale_x_discrete(labels = large_labels) +
  labs(x = NULL, y = 'Expression', caption = 'Letter-value plot') +
  theme_classic(base_size = 10) + theme(legend.position = 'none')

# 4. SPLIT VIOLIN: only complete clusters with >=30 observations in each ordered condition survive.
split_group <- interaction(df_paired$cluster, df_paired$condition, lex.order = TRUE)
split_bw <- safe_violin_bw(transform(df_paired, split_group = split_group), 'split_group', 'expression')
split_labels <- sprintf('%s\n(Control n=%d, Treatment n=%d)', levels(df_paired$cluster),
  as.integer(tapply(df_paired$condition == 'Control', df_paired$cluster, sum)),
  as.integer(tapply(df_paired$condition == 'Treatment', df_paired$cluster, sum)))
p_split <- ggplot(df_paired, aes(cluster, expression, fill = condition,
                                  group = interaction(cluster, condition, lex.order = TRUE))) +
  geom_split_violin(alpha = 0.7, trim = TRUE, bw = split_bw) +
  geom_boxplot(width = 0.15, position = position_dodge(0.5), outlier.shape = NA) +
  scale_fill_manual(values = c(Control = '#56B4E9', Treatment = '#D55E00')) +
  scale_x_discrete(labels = split_labels) +
  labs(x = 'Cluster', y = 'Expression', fill = NULL,
       caption = 'Incomplete or <30-per-condition clusters omitted') +
  theme_classic(base_size = 10) + theme(legend.position = 'top')

# 5. EXPORT: 89 mm is Nature single-column width.
ggsave('raincloud.pdf', p_raincloud, width = 89, height = 70, units = 'mm', device = cairo_pdf)
