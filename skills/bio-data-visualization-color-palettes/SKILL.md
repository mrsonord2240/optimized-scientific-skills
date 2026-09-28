---
name: bio-data-visualization-color-palettes
description: Select colormaps and qualitative palettes for scientific figures using perceptual-uniformity, color-vision-deficiency safety, and luminance-monotonicity criteria. Covers Crameri scientific colormaps, viridis/cividis/magma, Okabe-Ito categorical, ColorBrewer, and the rainbow/jet critique. Use when choosing palettes for heatmaps, scatter, networks, or any encoding where color carries quantitative or categorical meaning.
tool_type: mixed
primary_tool: viridis
goal_approach_exempt: true
---

## Version Compatibility

Examples were checked with viridis 0.6.5, RColorBrewer 1.1.3, scico 1.5.0, ggsci 5.2.0, colorspace 2.1.2, matplotlib 3.11.2, cmcrameri 1.10, and colorspacious 1.1.2.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Color Palettes for Scientific Visualization

**"Pick a color palette"** -> Choose a colormap that (a) is perceptually uniform along the relevant data axis, (b) remains interpretable under common color-vision deficiencies, (c) prints correctly to grayscale, and (d) matches the data type: sequential, diverging, cyclic, or qualitative.

- R: `viridis::viridis`, `scico::scale_color_scico`, `RColorBrewer::brewer.pal`
- Python: `matplotlib.colormaps`, `seaborn.color_palette`, `cmcrameri.cm`

## The Three Modern Standards

1. **Perceptual uniformity** -- equal data steps produce equal perceived color steps. viridis (van der Walt 2015), cividis (Nuñez 2018), and the Crameri family (batlow, roma, vik) are designed for this. Jet, rainbow, and red->green are not.

2. **Color vision deficiency (CVD) safety** -- ~6% of males have deuteranopia or protanopia (red-green deficiency). cividis was explicitly designed for CVD viewing (Nuñez 2018, *PLOS ONE* 13:e0199239). The Okabe-Ito qualitative palette (popularized in Wong 2011, *Nat Methods* 8:441) is the CVD-safe categorical default.

3. **Grayscale monotonicity** -- a perceptually uniform sequential colormap has monotonically changing luminance. Desaturate the actual palette; if its order remains readable, it is luminance-monotonic. This is not a substitute for CVD simulation.

## Palette Type by Data Type

| Data type | Use | Avoid |
|-----------|-----|-------|
| Sequential (expression, coverage, density) | viridis, magma, cividis, batlow, lipari | jet, rainbow, hsv |
| Diverging (log fold change, z-score, signed correlation) | vik, roma, RdBu, BrBG, PiYG; use symmetric limits | jet, rainbow |
| Cyclic (phase, time-of-day, angle) | romaO, vikO, twilight | linear sequential (wrap creates an artifactual jump) |
| Categorical (≤8 groups) | Named Okabe-Ito mapping | unnamed vectors; Set1 if CVD matters |
| Categorical (9-20 groups) | No single palette is CVD-safe; pair color with shape, direct labels, or facets | relying on hue alone |
| Categorical (>20) | Reconsider the design | adding more colors |

## The Crameri Scientific Colormaps

Crameri 2020 *Nat Commun* 11:5444 documented the prevalence of misleading palettes and released perceptually uniform scientific colormaps. Centres below were measured from 255 colors with scico 1.5.0; exact hexes can vary by implementation and sample count.

| Crameri name | Type | Use case or measured centre |
|--------------|------|-----------------------------|
| `batlow` | sequential | Default jet replacement; dark blue -> ochre -> light pink |
| `lipari` | sequential | Higher-saturation sequential alternative |
| `vik` | diverging | Blue -> warm off-white (`#EBE5E0`) -> red |
| `roma` | diverging | Red -> pale green (`#C0E9C2`) -> blue |
| `bam` | diverging | Magenta -> near-white (`#F5F0F0`) -> green |
| `romaO` | cyclic | Phase, time-of-day, angle data |
| `vikO` | cyclic | Diverging cyclic data |

```r
library(scico)
ggplot(df, aes(x, y, fill = value)) + geom_tile() +
    scale_fill_scico(palette = 'batlow')

vmax <- quantile(abs(df$lfc), 0.99, na.rm = TRUE)
ggplot(df, aes(x, y, fill = lfc)) + geom_tile() +
    scale_fill_scico(palette = 'vik', midpoint = 0,
                     limits = c(-vmax, vmax), oob = scales::squish)
```

```python
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
plt.imshow(data, cmap=cm.batlow)
vmax = np.nanpercentile(np.abs(data), 99)
plt.imshow(data, cmap=cm.vik, vmin=-vmax, vmax=vmax)
```

## viridis Family (matplotlib default since 2.0)

```r
library(viridis)
scale_color_viridis_c(option = 'viridis')   # dark blue -> yellow
scale_color_viridis_c(option = 'magma')     # black -> red -> yellow
scale_color_viridis_c(option = 'inferno')   # black -> purple -> yellow
scale_color_viridis_c(option = 'plasma')    # purple -> pink -> yellow
scale_color_viridis_c(option = 'cividis')   # CVD-optimized
```

```python
plt.imshow(data, cmap='viridis')   # also magma, inferno, plasma, cividis
```

Use viridis, cividis, or batlow as the default replacement for jet. Turbo is a smoother jet-like palette and can be a last resort when that appearance is required, but its luminance rises then falls, so it is neither perceptually uniform nor grayscale-monotonic.

## Okabe-Ito Categorical Palette (Wong 2011)

Name every mapping so subsetting or reordering factor levels cannot silently recolor groups. Reserve light grey separately rather than consuming one of the seven chromatic colors.

```r
cell_levels <- c('T', 'B', 'NK', 'Monocyte', 'Dendritic', 'Platelet', 'Erythroid')
okabe_ito <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442',
               '#0072B2', '#D55E00', '#CC79A7')
cell_colors <- c(setNames(okabe_ito, cell_levels), Unassigned = '#BBBBBB')
scale_color_manual(values = cell_colors, drop = FALSE)
```

Named vectors keep colors stable across panels and subsets. Under deutan simulation, `#CC79A7` and `#009E73` can approach mid-grey; check the full mapping and use a marker shape or direct label when the reserve class still collides. R 4.0+ also provides `palette.colors(8, 'Okabe-Ito')`. Matplotlib has no `colorblind` style; `seaborn-v0_8-colorblind` is only a six-color cycle, so use the explicit named mapping when identity must be stable.

For DE plots, the canonical assignment is Up = `#D55E00` (vermilion), Down = `#0072B2` (blue), NS = `#999999` (grey).

## ColorBrewer (Harrower & Brewer 2003)

| Palette | Type | Measured centre at maximum odd size |
|---------|------|-------------------------------------|
| `YlOrRd` | sequential | n/a |
| `RdBu` | diverging | near-neutral `#F7F7F7` at 6/11 |
| `Dark2` | qualitative | n/a; simulate before accessibility use |

```r
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 9, name = 'YlOrRd')
brewer.pal(n = 11, name = 'RdBu')
```

ColorBrewer's sequential and diverging palettes remain useful publication defaults. A near-neutral built-in centre is intentional; zero is anchored by `midpoint = 0` or symmetric limits, not by requiring the centre hex to be pure white.

## Scientific Journal Brand Palettes

```r
library(ggsci)
scale_color_npg()       # Nature Publishing Group
scale_color_aaas()      # Science (AAAS)
scale_color_lancet()    # Lancet
scale_color_jama()      # JAMA
scale_color_jco()       # JCO
scale_color_nejm()      # NEJM
```

These are CVD-imperfect. Use journal palettes for stylistic compliance, not as accessibility defaults, and simulate the chosen colors.

## CVD Simulation -- The Mandatory Check

`colorspace::cvd_emulator()` takes an image filename and launches an interactive emulator; it does not transform a palette. Use `deutan()`, `protan()`, and `tritan()` for palette vectors.

```r
library(colorspace)
simulated <- list(deutan = deutan(palette), protan = protan(palette),
                  tritan = tritan(palette))
demoplot(simulated$deutan, type = 'heatmap')
sapply(simulated, function(p) min(dist(coords(as(hex2RGB(p), 'LAB')))))
```

```python
import numpy as np
from colorspacious import cspace_convert
from matplotlib.colors import to_rgb
from scipy.spatial.distance import pdist
rgb = np.array([to_rgb(color) for color in palette])
for name in ('deuteranomaly', 'protanomaly', 'tritanomaly'):
    spec = {'name': 'sRGB1+CVD', 'cvd_type': name, 'severity': 100}
    sim = np.clip(cspace_convert(rgb, spec, 'sRGB1'), 0, 1)
    print(name, pdist(cspace_convert(sim, 'sRGB1', 'CAM02-UCS')).min())
```

Report the minimum pairwise distance and inspect the simulated figure. If categories become hard to distinguish, change the encoding; a single distance cutoff is not a guarantee of accessibility.

## Grayscale Monotonicity Test

```r
library(scales)
pal <- viridis::viridis(10)
show_col(pal)
show_col(desaturate(pal))
L <- coords(as(hex2RGB(pal), 'LAB'))[, 'L']
all(diff(L) >= 0) || all(diff(L) <= 0)
```

This tests the actual palette rather than replacing it with a generic 0-100 grey ramp. Grayscale monotonicity catches false bands in sequential maps; CVD simulation is a separate check. Rainbow, jet, and turbo fail the grayscale test because their luminance is non-monotonic.

## Diverging Palette Setup (LFC, z-score)

Built-in palettes such as vik, roma, bam, and RdBu use near-neutral, sometimes tinted centres by design. Anchor zero with `midpoint = 0` and symmetric limits. Reserve pure white for a custom ramp whose specification requires it.

```r
library(circlize)
col_fun <- colorRamp2(c(-2, 0, 2), c('#0072B2', 'white', '#D55E00'))
```

```python
import numpy as np
import matplotlib.pyplot as plt
vmax = np.nanpercentile(np.abs(data), 99)
plt.imshow(data, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
```

Do not use `vmin=data.min(), vmax=data.max()` for signed data: asymmetric bounds move zero away from the palette centre.

## Nonzero References and Missing Values

The meaningful centre of a diverging scale is the scientific reference, not necessarily zero. State that reference explicitly, choose bounds on both sides of it, and keep missing or non-finite cells out of the quantitative scale. Give missing values a separate neutral color and label it as missing rather than implying a numeric value.

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

reference = 0.5
masked = np.ma.masked_invalid(data)
norm = TwoSlopeNorm(vmin=0.0, vcenter=reference, vmax=1.0)
cmap = plt.get_cmap('RdBu_r').copy()
cmap.set_bad('#BBBBBB')
plt.imshow(masked, cmap=cmap, norm=norm)
```

`TwoSlopeNorm` maps the declared reference to the palette centre. Validate that `vmin < reference < vmax`; if that is not true, a diverging scale is not defined for the supplied bounds. Report how many values were masked, and show missingness in the legend or caption.

## Palette Audit Response Contract

When recommending or auditing a palette, report this compact checklist:

- palette name and type: sequential, diverging, cyclic, or categorical
- normalization, bounds, and reference: include clipping or percentile rules
- missing values: masked count and their separate neutral encoding
- CVD minimum distances: deutan, protan, and tritan for categorical mappings
- luminance verdict: monotonic, non-monotonic, or not applicable
- redundant encoding: shape, labels, facets, or none with a reason
- visual inspection: artifact opened, simulation inspected, and any limitation

## Custom Palette Construction

```r
my_palette <- c(Control = '#0072B2', Treatment = '#D55E00', Vehicle = '#009E73')
scale_color_manual(values = my_palette)

# Use an odd count when an exact sampled white centre is required.
colorRampPalette(c('#0072B2', 'white', '#D55E00'))(101)
```

```python
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('cvd_div', ['#0072B2', '#FFFFFF', '#D55E00'])
```

For a custom diverging palette, choose endpoints with similar luminance so neither side dominates, specify its centre deliberately, and run both grayscale and CVD checks.

## Common Failure Modes

### Categorical palette with too many colors

Human color discrimination saturates around 8-10 hues, and the tested 9-20 color palettes had close pairs under deutan simulation. Facet, aggregate small groups into "Other," add marker shape, or use direct labels.

### Rainbow / jet still in use

Rainbow has non-monotonic luminance and creates artifactual yellow bands. Migrate sequential data to viridis or batlow and diverging data to vik or roma. Turbo preserves the jet-like appearance but also fails grayscale monotonicity.

### Wrong cmap for non-numeric data

A continuous gradient falsely implies order for nominal cluster IDs. Use a named qualitative mapping, preserving identities across panels and subsets.

### CVD-unsafe palette in a CVD-sensitive figure

Pre-flight the final mapping with the CVD simulation above. Switch to Okabe-Ito for up to eight categories or cividis for sequential data, and add a redundant encoding when colors still collide.

## References

- Crameri F, Shephard GE, Heron PJ. 2020. The misuse of colour in science communication. *Nat Commun* 11:5444. doi:10.1038/s41467-020-19160-7
- Harrower M, Brewer CA. 2003. ColorBrewer.org: an online tool for selecting colour schemes for maps. *Cartogr J* 40(1):27-37.
- Nuñez JR, Anderton CR, Renslow RS. 2018. Optimizing colormaps with consideration for color vision deficiency to enable accurate interpretation of scientific data. *PLOS ONE* 13(7):e0199239.
- Borland D, Taylor RM II. 2007. Rainbow color map (still) considered harmful. *IEEE Comput Graph Appl* 27(2):14-17.
- Light A, Bartlein PJ. 2004. The end of the rainbow? Color schemes for improved data graphics. *Eos Trans AGU* 85(40):385,391.
- Wong B. 2010. Points of view: Color coding. *Nat Methods* 7(8):573.
- Wong B. 2011. Points of view: Color blindness. *Nat Methods* 8(6):441.
- Gehlenborg N, Wong B. 2012. Points of view: Mapping quantitative data to color. *Nat Methods* 9(8):769.

## Related Skills

- data-visualization/heatmaps-clustering - Robust diverging bounds for heatmaps
- data-visualization/volcano-and-ma-plots - Okabe-Ito Up/Down/NS conventions
- data-visualization/ggplot2-fundamentals - Applying palettes in ggplot2 scales
- data-visualization/dimensionality-reduction-plots - Categorical palettes for cluster labels
