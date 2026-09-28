# Purpose: self-contained examples of sequential, diverging, and stable categorical palettes.
# Inputs: simulated PCA-like data generated below. Usage: Rscript palette_examples.R
# Reference: ggplot2 4.0.3, viridis 0.6.5, scico 1.5.0, patchwork 1.3.2.
library(ggplot2)
library(viridis)
library(scico)
library(patchwork)

set.seed(42)
n_samples <- 100
group_levels <- c('Control', 'DrugA', 'DrugB', 'DrugC', 'DrugD')
group <- factor(sample(group_levels, n_samples, replace = TRUE), levels = group_levels)
group_centers <- list(Control = c(0, 0), DrugA = c(3, 1), DrugB = c(-2, 2),
                      DrugC = c(1, -3), DrugD = c(-3, -2))
df <- data.frame(
    PC1 = vapply(group, function(g) rnorm(1, group_centers[[g]][1], 0.8), numeric(1)),
    PC2 = vapply(group, function(g) rnorm(1, group_centers[[g]][2], 0.8), numeric(1)),
    expression = rexp(n_samples, rate = 0.5),
    log2FC = rnorm(n_samples, 0, 1.5),
    group = group
)

vmax <- quantile(abs(df$log2FC), 0.99, na.rm = TRUE)
okabe_ito <- setNames(
    c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2'),
    group_levels
)

p1 <- ggplot(df, aes(PC1, PC2, color = expression)) +
    geom_point(size = 3) +
    scale_color_scico(palette = 'batlow') +
    labs(title = 'Batlow (Sequential)', color = 'Expression') +
    theme_minimal()

p2 <- ggplot(df, aes(PC1, PC2, fill = log2FC)) +
    geom_point(shape = 21, color = 'grey20', stroke = 0.25, size = 3) +
    scale_fill_scico(palette = 'vik', midpoint = 0,
                     limits = c(-vmax, vmax), oob = scales::squish) +
    labs(title = 'Vik (Diverging, Symmetric)', fill = 'log2FC') +
    theme_minimal() +
    theme(panel.background = element_rect(fill = 'grey95', color = NA))

p3 <- ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = okabe_ito, drop = FALSE) +
    labs(title = 'Named Okabe-Ito', color = 'Treatment') +
    theme_minimal()

p4 <- ggplot(subset(df, group != 'Control'), aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = okabe_ito, drop = FALSE) +
    labs(title = 'Same Mapping After Subsetting', color = 'Treatment') +
    theme_minimal()

combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = 'Color Palette Examples',
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 14,
                                                            face = 'bold')))

ggsave('palette_examples.pdf', combined, width = 12, height = 10)
message('Saved: palette_examples.pdf')
cat('\nRecommended R palettes:\n')
cat('Sequential: viridis, magma, cividis, batlow, lipari\n')
cat('Diverging: vik, roma, RdBu\n')
cat('Categorical: named Okabe-Ito mapping for up to 8 groups\n')
