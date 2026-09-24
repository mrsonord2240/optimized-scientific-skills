## SILAC Quantification

**Goal:** Compute heavy/light ratios while preserving on/off biology and flagging label artifacts.

**Approach:** A protein present only in the heavy channel may be the interesting biology (or a detection-limit dropout), so do not discard it -- but keep it in a presence flag, not as +/-Inf inside the ratio matrix, where it turns pandas SDs into NaN and stops `eBayes(trend=TRUE, robust=TRUE)`. Verify labeling efficiency (>=95%, target 97-98%) on a heavy-only pilot and assess Arg->Pro conversion before trusting any ratio -- both bias every ratio in the same direction, so neither shows up as extra scatter.

### Check labeling efficiency and Arg->Pro first
```bash
python scripts/silac_checks.py pilot_evidence.txt   # heavy-only pilot peptide table (tab-separated) -> JSON
```
Or `from silac_checks import silac_labeling_efficiency, arg_to_pro_shift`. When the usual `Sequence` column is present, `silac_labeling_efficiency` excludes Pro-containing peptides before reporting intensity-weighted incorporation: Arg->Pro conversion otherwise drains their heavy signal and makes incorporation read falsely low. Its JSON reports how many signal peptides it used and excluded; without `Sequence`, `sequence_column_used` is false, so use the all-peptide estimate only after separately ruling out conversion. It also reports peptides below 95% and the log2 H/L bias a 1:1 mix would show (`eff / (2 - eff)`, -0.20 at 93%: the unincorporated fraction of the heavy sample is counted in the LIGHT channel). `arg_to_pro_shift` regresses log2 H/L on proline count: the dose slope is what makes it Arg->Pro (a flat offset is incomplete labeling); the direct route is to re-search the pilot with Pro6 variable and take I(Pro6)/(I(Pro6)+I(Pro0)).

Seeded pilot recovering the planted 0.93 incorporation and 0.08 conversion: `examples/lfq_normalization.py`.

### Ratios that keep on/off biology
```python
import numpy as np

# Arg10/Lys8 is the common pairing (avoids overlap with the +6 isotope envelope)
SILAC_SHIFTS = {'Arg10': 10.008269, 'Lys8': 8.014199, 'Arg6': 6.020129, 'Lys6': 6.020129}

def silac_log2_ratio(heavy, light):
    '''Vectorized over arrays/Series: log2 H/L (NaN unless both channels quantified) plus a presence flag.'''
    heavy, light = np.asarray(heavy, dtype=float), np.asarray(light, dtype=float)
    h, l = heavy > 0, light > 0    # NaN compares False
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(h & l, np.log2(heavy / light), np.nan)
    presence = np.select([h & l, h, l], ['both', 'H-only', 'L-only'], default='none')    # report H-only/L-only separately
    return ratio, presence
```
