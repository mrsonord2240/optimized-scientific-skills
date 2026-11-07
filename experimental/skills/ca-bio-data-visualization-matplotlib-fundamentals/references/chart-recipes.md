# matplotlib chart and axis recipes

Tested on matplotlib 3.11.2. `Axes.boxplot(labels=)` was renamed `tick_labels=` (3.9) and is rejected on 3.11.

## Common chart types

```python
# Scatter -- always rasterized for >1000 points
ax.scatter(x, y, c=values, cmap='viridis', s=8, alpha=0.6,
           edgecolors='none', rasterized=True)
fig.colorbar(ax.collections[0], ax=ax, label='Expression', shrink=0.8)

# Line
ax.plot(x, y1, color='#0072B2', label='Control', linewidth=1)
ax.plot(x, y2, color='#D55E00', label='Treatment', linewidth=1)
ax.fill_between(x, y_low, y_high, color='#0072B2', alpha=0.2)
ax.legend(frameon=False, fontsize=6)

# Bar
ax.bar(categories, values, color='#0072B2', edgecolor='black', linewidth=0.5)

# Box / violin (prefer seaborn for these -- see distribution-plots)
ax.boxplot([group_a, group_b, group_c], tick_labels=['A', 'B', 'C'],
           patch_artist=True, boxprops=dict(facecolor='#0072B2', alpha=0.7))

# Histogram
ax.hist(values, bins=30, color='#0072B2', edgecolor='white', linewidth=0.5)

# Heatmap (prefer seaborn for clustered; see heatmaps-clustering)
im = ax.imshow(matrix, cmap='RdBu_r', aspect='auto', vmin=-vmax, vmax=vmax)
fig.colorbar(im, ax=ax, label='Z-score')
```

## Axis formatting

```python
# Log scale
ax.set_yscale('log')

# Scientific notation
from matplotlib.ticker import ScalarFormatter
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))

# Date axis
import matplotlib.dates as mdates
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

# Tick frequency
import numpy as np
ax.set_xticks(np.arange(0, 6, 2))                  # one tick per label: 0, 2, 4
ax.set_xticklabels(['A', 'B', 'C'], rotation=45, ha='right')

# Grid
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
```
