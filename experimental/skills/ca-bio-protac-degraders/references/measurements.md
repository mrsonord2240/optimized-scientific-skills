## Cooperativity (Alpha)

Cooperativity quantifies how the ternary complex stabilizes (or destabilizes) the binary binding:

```
alpha = (Kd_binary,target) / (Kd_ternary,target)
```

- alpha > 1: positive cooperativity (ternary stronger than binary)
- alpha = 1: no cooperativity (independent binding)
- alpha < 1: negative cooperativity (mutual destabilization)

Positive cooperativity can favor ternary-complex formation, but the preferred alpha is system- and assay-dependent and does not alone establish degradation efficacy. Alpha must be measured from binary and ternary binding experiments; a predicted structure does not directly provide it.

Measure with ITC (isothermal titration calorimetry) or SPR/BLI titrations of binary vs ternary.

Runnable: `scripts/cooperativity_dc50.py`'s `cooperativity_alpha(kd_binary, kd_ternary)` applies
this formula to your measured Kd pair and labels positive/negative/no cooperativity.

## DC50 / Dmax Characterization

In cellular assays:
- **DC50**: PROTAC concentration for 50% degradation (analogous to IC50)
- **Dmax**: maximum fraction degraded at any concentration

| Property | What to report | Interpretation |
|----------|----------------|----------------|
| DC50 | Concentration producing 50% of the assay's fitted maximal degradation | Compare only across matched assay conditions; no universal clinical cutoff |
| Dmax | Maximum observed or fitted degradation and uncertainty | Required depletion is target- and phenotype-dependent |
| Hook effect | Full concentration-response range and concentration of any downturn | A high-concentration effect; its location is system- and assay-dependent |
| Cooperativity | Alpha from matched binary and ternary binding experiments | Distinct from the hook effect and insufficient by itself to predict degradation |

**Hook effect**: at high PROTAC concentrations, binary complexes (PROTAC-target alone, PROTAC-E3 alone) dominate, and ternary complex formation drops. Dose-response curves are bell-shaped.

Runnable: `scripts/cooperativity_dc50.py` fits DC50/Dmax with a 3-parameter Hill curve, flags a
hook effect from a >15-percentage-point downturn after the peak, and -- when a hook is flagged --
restricts the DC50/Dmax fit to the ascending arm rather than fitting a monotonic sigmoid to
non-monotonic data. When the peak sits under 10x the fitted DC50 (hook onset close to DC50), the
ascending arm never reaches its plateau and Dmax is underestimated (15 pp on a noise-free
synthetic curve); `fit_dc50()` then returns `plateau_reached=False` and a `caveat` -- report Dmax as
a lower bound. It ships a seeded synthetic dose-response curve and checks the fit recovers
the curve's own planted DC50/Dmax within tolerance, and that a narrow-separation hook is caveated.
