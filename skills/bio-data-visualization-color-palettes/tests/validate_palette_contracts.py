"""Regression-test Python palette guidance on the audit spatial array.

Usage: py.sh validate_palette_contracts.py <audit-data-dir>
"""

from pathlib import Path
import py_compile
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np
from colorspacious import cspace_convert
from matplotlib.colors import TwoSlopeNorm, to_rgb
from scipy.spatial.distance import pdist


if len(sys.argv) != 2:
    raise SystemExit("usage: validate_palette_contracts.py <audit-data-dir>")
py_compile.compile(__file__, cfile=str(Path(tempfile.gettempdir()) / "validate_palette_contracts.pyc"),
                   doraise=True)
data_dir = Path(sys.argv[1])
skill_text = (Path(__file__).parent.parent / "SKILL.md").read_text(encoding="utf-8")
assert "## Nonzero References and Missing Values" in skill_text
assert "## Palette Audit Response Contract" in skill_text
for required_field in (
    "palette name and type",
    "normalization, bounds, and reference",
    "CVD minimum distances",
    "luminance verdict",
    "redundant encoding",
    "visual inspection",
):
    assert required_field in skill_text

data = np.load(data_dir / "synthetic_spatial_expr.npy")
vmax = np.nanpercentile(np.abs(data), 99)
assert np.isfinite(vmax) and vmax > 0

figure, axis = plt.subplots()
image = axis.imshow(data, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
assert image.norm.vmin == -image.norm.vmax
plt.close(figure)

reference = 0.5
reference_data = np.array([[0.0, reference, 1.0], [np.nan, np.inf, -np.inf]])
masked = np.ma.masked_invalid(reference_data)
norm = TwoSlopeNorm(vmin=0.0, vcenter=reference, vmax=1.0)
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("#BBBBBB")
assert norm(reference) == 0.5
assert int(masked.mask.sum()) == 3
assert np.allclose(cmap.get_bad(), np.array([187, 187, 187, 255]) / 255)

palette = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#BBBBBB",
]
rgb = np.array([to_rgb(color) for color in palette])
minimum_distances = {}
for name in ("deuteranomaly", "protanomaly", "tritanomaly"):
    spec = {"name": "sRGB1+CVD", "cvd_type": name, "severity": 100}
    simulated = np.clip(cspace_convert(rgb, spec, "sRGB1"), 0, 1)
    distance = pdist(cspace_convert(simulated, "sRGB1", "CAM02-UCS")).min()
    assert np.isfinite(distance) and distance > 0
    minimum_distances[name] = float(distance)

samples = np.linspace(0, 1, 256)
viridis_rgb = plt.get_cmap("viridis")(samples)[:, :3]
turbo_rgb = plt.get_cmap("turbo")(samples)[:, :3]
viridis_l = cspace_convert(viridis_rgb, "sRGB1", "CIELab")[:, 0]
turbo_l = cspace_convert(turbo_rgb, "sRGB1", "CIELab")[:, 0]
assert np.all(np.diff(viridis_l) >= 0)
assert not (np.all(np.diff(turbo_l) >= 0) or np.all(np.diff(turbo_l) <= 0))

print("Python palette contracts PASS")
print("shape:", data.shape, "q99:", vmax)
print("CVD minimum CAM02-UCS distances:", minimum_distances)
