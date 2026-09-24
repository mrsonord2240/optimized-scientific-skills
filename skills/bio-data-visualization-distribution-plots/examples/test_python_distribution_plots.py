# Purpose: smoke-test the current Python distribution-plot snippets.
# Inputs: none. Usage: python test_python_distribution_plots.py

from pathlib import Path
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ptitprince as pt
import seaborn as sns

rng = np.random.default_rng(20260923)
frame = pd.DataFrame({
    "group": np.repeat(["Control", "Treated"], 40),
    "value": np.r_[rng.normal(5, 0.7, 40), rng.normal(6, 0.8, 40)],
})
palette = {"Control": "#0072B2", "Treated": "#D55E00"}
out_dir = Path(tempfile.mkdtemp(prefix="distribution-plots-python-"))

fig, ax = plt.subplots(figsize=(5, 3.5))
pt.RainCloud(x="group", y="value", data=frame, hue="group", palette=palette,
             bw="scott", cut=0, width_viol=0.6, orient="h", ax=ax)
raincloud_path = out_dir / "raincloud.png"
fig.savefig(raincloud_path, dpi=120)
plt.close(fig)

fig, ax = plt.subplots(figsize=(4, 3.5))
sns.boxenplot(data=frame, x="group", y="value", hue="group", palette=palette,
              legend=False, ax=ax)
boxen_path = out_dir / "boxen.png"
fig.savefig(boxen_path, dpi=120)
plt.close(fig)

assert raincloud_path.stat().st_size > 2000
assert boxen_path.stat().st_size > 2000
print(f"PASS: {raincloud_path.name}={raincloud_path.stat().st_size}; {boxen_path.name}={boxen_path.stat().st_size}")
