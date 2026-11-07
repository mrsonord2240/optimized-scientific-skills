"""Plot two CSV columns as an 89 x 70 mm publication scatter."""
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input_csv')
parser.add_argument('x')
parser.add_argument('y')
parser.add_argument('output')
args = parser.parse_args()
df = pd.read_csv(args.input_csv)
mpl.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 7, 'savefig.dpi': 300})
fig, ax = plt.subplots(figsize=(89 / 25.4, 70 / 25.4), layout='constrained')
ax.scatter(df[args.x], df[args.y], c='#0072B2', s=10, alpha=0.7, edgecolors='none', rasterized=True)
ax.set_xlabel(args.x)
ax.set_ylabel(args.y)
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(args.output)
plt.close(fig)
