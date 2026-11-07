# Purpose: self-contained palette selection examples for four scientific data types.
# Inputs: simulated tile, point, and differential-expression data generated below.
# Usage: Rscript palettes_phd.R
# Reference: scico 1.5.0, viridis 0.6.5, colorspace 2.1.2, ggplot2 4.0.3.

library(ggplot2)
library(scico)
library(viridis)
library(colorspace)

set.seed(42)
df <- expand.grid(x = 1:20, y = 1:20)
df$expression <- rexp(nrow(df), rate = 0.5)
df$lfc <- rnorm(nrow(df), 0, 1.5)
df$phase <- atan2(df$y - mean(df$y), df$x - mean(df$x)) %% (2 * pi)
cell_levels <- c('T', 'B', 'NK', 'Monocyte', 'Dendritic', 'Platelet',
                 'Erythroid', 'Unassigned')
df$cell_type <- factor(sample(cell_levels, nrow(df), replace = TRUE),
                       levels = cell_levels)
df$group <- factor(sample(LETTERS[1:5], nrow(df), replace = TRUE))

de_df <- data.frame(log2FC = rnorm(500, 0, 1.5), pvalue = runif(500))
de_df$neg_log10_p <- -log10(de_df$pvalue)
de_df$significance <- factor(
    ifelse(de_df$log2FC > 1 & de_df$pvalue < 0.05, 'Up',
           ifelse(de_df$log2FC < -1 & de_df$pvalue < 0.05, 'Down', 'NS')),
    levels = c('Up', 'Down', 'NS')
)

# 1. SEQUENTIAL -- Crameri batlow or viridis cividis
ggplot(df, aes(x, y, fill = expression)) + geom_tile() +
    scale_fill_scico(palette = 'batlow')

# 2. DIVERGING -- Crameri vik with symmetric robust bounds
vmax <- quantile(abs(df$lfc), 0.99, na.rm = TRUE)
ggplot(df, aes(x, y, fill = lfc)) + geom_tile() +
    scale_fill_scico(palette = 'vik', midpoint = 0,
                     limits = c(-vmax, vmax), oob = scales::squish)

# 3. CYCLIC -- phase, time-of-day, or angle
ggplot(df, aes(x, y, color = phase)) + geom_point() +
    scale_color_scico(palette = 'romaO')

# 4. CATEGORICAL <=8 -- named Okabe-Ito mapping plus reserved light grey
okabe_ito <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442',
               '#0072B2', '#D55E00', '#CC79A7')
cell_colors <- c(setNames(okabe_ito, cell_levels[1:7]), Unassigned = '#BBBBBB')
ggplot(df, aes(x, y, color = cell_type)) + geom_point() +
    scale_color_manual(values = cell_colors, drop = FALSE)

# 5. CATEGORICAL DE convention -- Up/Down/NS
de_colors <- c(Up = '#D55E00', Down = '#0072B2', NS = '#999999')
ggplot(de_df, aes(log2FC, neg_log10_p, color = significance)) + geom_point() +
    scale_color_manual(values = de_colors)

# 6. CVD SIMULATION -- mandatory and separate from grayscale testing
palette <- scico(8, palette = 'batlow')
demoplot(palette, type = 'heatmap')
demoplot(deutan(palette), type = 'heatmap')
demoplot(protan(palette), type = 'heatmap')
demoplot(tritan(palette), type = 'heatmap')
sapply(list(deutan = deutan(palette), protan = protan(palette),
            tritan = tritan(palette)),
       function(p) min(dist(coords(as(hex2RGB(p), 'LAB')))))

# 7. GRAYSCALE MONOTONICITY TEST
library(scales)
show_col(palette)
show_col(desaturate(palette))
L <- coords(as(hex2RGB(palette), 'LAB'))[, 'L']
stopifnot(all(diff(L) >= 0) || all(diff(L) <= 0))

# 8. COMPARE -- the jet/rainbow trap visualized
show_col(rainbow(10))
show_col(desaturate(rainbow(10)))
show_col(viridis(10))
show_col(desaturate(viridis(10)))

# 9. JOURNAL-BRAND PALETTES -- stylistic, not accessibility defaults
library(ggsci)
ggplot(df, aes(x, y, color = group)) + geom_point() + scale_color_npg()
