"""Render an MA plot from a CSV of finite shrunken effects and positive means."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ma_plot import ma_plot

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input_csv')
parser.add_argument('output')
args = parser.parse_args()
df = pd.read_csv(args.input_csv)
if not {'baseMean', 'log2FoldChange', 'padj'} <= set(df.columns):
    raise ValueError('missing required columns')
if not np.isfinite(df['baseMean']).all() or (df['baseMean'] <= 0).any():
    raise ValueError('baseMean must be finite and positive')
if not np.isfinite(df['log2FoldChange']).all():
    raise ValueError('effects must be finite')
p = df['padj'].dropna()
if not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
    raise ValueError('nonmissing padj must be in [0, 1]')
ax = ma_plot(df)
ax.figure.savefig(args.output, dpi=300)
plt.close(ax.figure)
