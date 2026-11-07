args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop('usage: plot_palette.R data.csv x y value mode output.pdf')
library(ggplot2)
library(scico)
df <- read.csv(args[1], check.names = FALSE)
if (!all(args[2:4] %in% names(df))) stop('missing column')
names(df)[match(args[2:4], names(df))] <- c('.x', '.y', '.value')
mode <- args[5]
p <- ggplot(df, aes(x = .x, y = .y, fill = .value)) + geom_tile()
if (mode == 'sequential') {
    p <- p + scale_fill_scico(palette = 'batlow', na.value = '#BBBBBB')
} else if (mode == 'cyclic') {
    p <- p + scale_fill_scico(palette = 'romaO', na.value = '#BBBBBB')
} else if (mode == 'diverging') {
    v <- quantile(abs(df$.value), 0.99, na.rm = TRUE)
    if (!is.finite(v) || v <= 0) stop('nonzero finite range required')
    p <- p + scale_fill_scico(palette = 'vik', midpoint = 0, limits = c(-v, v),
                              oob = scales::squish, na.value = '#BBBBBB')
} else if (mode == 'categorical') {
    levels <- sort(unique(as.character(na.omit(df$.value))))
    if (length(levels) > 7 || length(levels) == 0) stop('requires one to seven categories')
    colors <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')
    df$.value <- factor(df$.value, levels = levels)
    p <- ggplot(df, aes(x = .x, y = .y, fill = .value)) + geom_tile() +
        scale_fill_manual(values = setNames(colors[seq_along(levels)], levels), na.value = '#BBBBBB', drop = FALSE)
} else stop('unknown mode')
ggsave(args[6], p, width = 8, height = 6)
