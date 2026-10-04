"""Reference: matplotlib 3.8+ (run clean on 3.11.2), seaborn 0.13+ | Verify API if version differs

PhD-level matplotlib pattern encoding the four correctness traps:
(1) Type-42 font embedding, (2) OO Figure/Axes API, (3) constrained_layout,
(4) rasterized point layer for vector axes + raster cells, (5) colours bound to category names.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

# Illustrative stand-ins so the recipes below run as-is; replace with real results in practice.
rng = np.random.default_rng(0)
df = pd.DataFrame({'pc1': rng.normal(size=500), 'pc2': rng.normal(size=500),
                   'cluster': rng.integers(0, 4, 500), 'log_fc': rng.normal(size=500),
                   'neg_log_p': rng.exponential(2, size=500),
                   'significance': rng.choice(['Up', 'Down', 'NS'], 500)})
panel_data = {k: {'x': rng.normal(size=200), 'y': rng.normal(size=200)} for k in list('abcdef')}
matrix = rng.normal(size=(30, 12))

# CVD-safe Okabe-Ito colours; categories are bound by name, never by data order.
okabe_ito = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']
sig_palette = {'NS': '#999999', 'Down': '#0072B2', 'Up': '#D55E00'}
sig_order = ['NS', 'Down', 'Up']
cluster_colors = [okabe_ito[i] for i in (0, 1, 2, 4)]   # skips yellow, faint on white

# 1. RCPARAMS FOR PUBLICATION COMPLIANCE
mpl.rcParams.update({
    'pdf.fonttype': 42,                          # TrueType -- searchable PDFs
    'ps.fonttype': 42,                           # TrueType EPS
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,                              # Nature 5-7 pt body
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'savefig.dpi': 300,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'lines.linewidth': 1.0,
})

# No savefig.bbox='tight': it crops to the artists and changes the page size away from the stated mm.
# 2. NATURE SINGLE-COLUMN SCATTER (89mm wide)
fig, ax = plt.subplots(figsize=(89 / 25.4, 70 / 25.4),
                        layout='constrained')
ax.scatter(df['pc1'], df['pc2'],
            c=df['cluster'].astype('category').cat.codes,
            cmap=ListedColormap(cluster_colors), s=8, alpha=0.7,
            edgecolors='none', rasterized=True)        # CRITICAL: raster for N>1000
ax.set_xlabel('PC1 (45.2%)')                            # variance-labeled
ax.set_ylabel('PC2 (12.1%)')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig('pca.pdf')                                  # Type-42 fonts; rasterized cells

# 3. MULTI-PANEL 2x3 GRID
fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 110 / 25.4),
                          layout='constrained')
panel_labels = list('abcdef')
for ax, label, (key, panel) in zip(axes.flat, panel_labels, panel_data.items()):
    ax.scatter(panel['x'], panel['y'], s=4, alpha=0.6, rasterized=True)
    ax.spines[['top', 'right']].set_visible(False)
    ax.text(-0.15, 1.05, label, transform=ax.transAxes,
            fontsize=7, fontweight='bold', va='top')    # panel label, kept inside the 5-7 pt range

fig.savefig('multipanel.pdf')

# 4. HEATMAP with Crameri batlow + raster
from cmcrameri import cm as cmc
fig, ax = plt.subplots(figsize=(89 / 25.4, 90 / 25.4),
                        layout='constrained')
vmax = np.quantile(np.abs(matrix), 0.99)               # robust 99th percentile bound
im = ax.imshow(matrix, cmap=cmc.vik, vmin=-vmax, vmax=vmax,
                aspect='auto', rasterized=True)
ax.set_xlabel('Sample')
ax.set_ylabel('Gene')
cbar = fig.colorbar(im, ax=ax, shrink=0.6, aspect=20, label='Z-score')
fig.savefig('heatmap.pdf')

# 5. SEABORN INTEGRATION -- shares matplotlib Figure/Axes
import seaborn as sns
fig, ax = plt.subplots(figsize=(89 / 25.4, 70 / 25.4),
                        layout='constrained')
sns.scatterplot(data=df, x='log_fc', y='neg_log_p',
                 hue='significance', hue_order=sig_order, palette=sig_palette,
                 s=8, alpha=0.7, ax=ax, rasterized=True,
                 legend='brief')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig('volcano_sns.pdf')

# 6. FONT EMBED VERIFICATION
# Run after savefig:
#   pdffonts pca.pdf
# Expected: "TrueType" or "Type 1" rows; NO "Type 3" rows. Page sizes equal the figsize in mm.

# 7. SEABORN OBJECTS GRAMMAR (ggplot-like; matplotlib 3.7+, seaborn 0.13+)
import seaborn.objects as so
# Plot's figure legend is placed outside the axes and falls off an exact-size page, so reserve the right 22% for it
# (extent) and move it inside the page after plot(). Limits: p._figure is a private seaborn attribute (verified on seaborn
# 0.13.2 only), and the fixed 22% fits short legend text only; a 29-character title covered 13 mm of the data. Use a smaller
# extent[2] or a shorter title and check the legend bbox, or the public route: Plot.on(fig), move fig.legends[0],
# fig.subplots_adjust(right=0.76) (measured inside the page).
plot = (so.Plot(df, x='log_fc', y='neg_log_p')
          .add(so.Dots(pointsize=2), color='significance')
          .scale(color=sig_palette)
          .theme({**sns.axes_style('ticks'), **mpl.rcParams})   # Plot's own theme otherwise overrides rcParams
          .layout(size=(89 / 25.4, 70 / 25.4), extent=[0, 0, 0.78, 1])
          .label(x='log2 fold change', y='-log10(p)'))
p = plot.plot()
legend = p._figure.legends[0]
legend.set_loc('center right')
legend.set_bbox_to_anchor((1.0, 0.5), transform=p._figure.transFigure)
p.save('volcano_so.pdf')
