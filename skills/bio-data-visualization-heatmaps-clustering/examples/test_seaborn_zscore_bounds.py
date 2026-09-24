"""Regression test: clustermap bounds must describe the row-z-scored matrix."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import seaborn as sns

df = pd.DataFrame(
    [[1.0, 2.0, 3.0, 7.0], [5.0, 7.0, 11.0, 13.0], [2.0, 2.5, 4.5, 8.0]],
    index=["gene1", "gene2", "gene3"],
    columns=["sample1", "sample2", "sample3", "sample4"],
)
row_sd = df.std(axis=1, ddof=1).replace(0, np.nan)
df_z = df.sub(df.mean(axis=1), axis=0).div(row_sd, axis=0).dropna(axis=0)
vmax = np.nanquantile(np.abs(df_z.to_numpy()), 0.99)

g = sns.clustermap(
    df_z,
    cmap="RdBu_r",
    center=0,
    vmin=-vmax,
    vmax=vmax,
    method="ward",
    metric="euclidean",
    z_score=None,
    figsize=(4, 4),
)
g.savefig("seaborn_zscore_bounds.png", dpi=150)

restored_data = g.data2d.reindex(index=df_z.index, columns=df_z.columns)
mesh_norm = g.ax_heatmap.collections[0].norm
assert np.allclose(restored_data.to_numpy(), df_z.to_numpy())
assert np.isclose(mesh_norm.vmin, -vmax)
assert np.isclose(mesh_norm.vmax, vmax)
assert np.isclose(vmax, np.nanquantile(np.abs(g.data2d.to_numpy()), 0.99))

print("PASS: seaborn renders the explicit row-z-scored data with matching robust bounds")
