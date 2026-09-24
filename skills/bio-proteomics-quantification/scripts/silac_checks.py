"""SILAC pilot checks: labeling efficiency and Arg->Pro conversion, on a HEAVY-ONLY pilot.

Input : tab-separated peptide table of the heavy-only pilot (MaxQuant evidence.txt or peptides.txt) with
        columns 'Intensity H', 'Intensity L' (efficiency) and 'Sequence', 'Ratio H/L' (Arg->Pro).
Output: JSON with both checks (Arg->Pro is skipped when its columns are absent).
Usage : python scripts/silac_checks.py pilot_evidence.txt
Import: from silac_checks import silac_labeling_efficiency, arg_to_pro_shift
Checked: numpy 2.5.3, pandas 3.0.5
"""
import numpy as np, pandas as pd

def silac_labeling_efficiency(pilot, heavy='Intensity H', light='Intensity L', sequence='Sequence'):
    '''Incorporation on a HEAVY-ONLY pilot (cells grown in heavy medium, NOTHING mixed in): every
    light ion there is unlabeled protein. pilot = the pilot's peptide table (MaxQuant evidence.txt).'''
    h = pd.to_numeric(pilot[heavy], errors='coerce').fillna(0.0)
    l = pd.to_numeric(pilot[light], errors='coerce').fillna(0.0)
    ok = (h + l) > 0
    if not ok.any():
        raise ValueError('no peptide has signal in either channel -- wrong columns or wrong file')
    # Arg->Pro conversion drains the heavy signal once per proline.  Including those peptides
    # makes a heavy-only pilot look incompletely labelled, so estimate incorporation on Pro-free
    # peptides when the usual MaxQuant Sequence column is available.  Keep the fallback for a
    # table without sequences, but make that limitation visible in the result.
    has_sequence = sequence in pilot.columns
    pro_free = ~pilot[sequence].fillna('').astype(str).str.contains('P', case=False, regex=False) if has_sequence else ok
    use = ok & pro_free
    if not use.any():
        raise ValueError('no Pro-free peptide has signal -- cannot separate labeling efficiency from Arg-to-Pro conversion')
    per_pep = h[use] / (h[use] + l[use])
    eff = float(h[use].sum() / (h[use] + l[use]).sum())    # intensity-weighted = the number to report
    # A 1:1 forward mix of these cells does NOT read log2 H/L = 0: the unincorporated (1 - eff) of
    # the heavy sample is counted in the LIGHT channel, so H/L = eff / (2 - eff) -- -0.20 at 93%.
    return {'n_peptides': int(use.sum()), 'n_signal_peptides': int(ok.sum()),
            'pro_containing_excluded': int((ok & ~pro_free).sum()) if has_sequence else 0,
            'sequence_column_used': has_sequence, 'incorporation': round(eff, 4),
            'median_peptide_incorporation': round(float(per_pep.median()), 4),
            'peptides_below_95pct': int((per_pep < 0.95).sum()),
            'expected_log2_HL_bias_at_1to1': round(float(np.log2(eff / (2 - eff))), 4),
            'pass_95pct': bool(eff >= 0.95)}

def arg_to_pro_shift(peptides, seq='Sequence', ratio='Ratio H/L'):
    '''Arg->Pro drains the heavy channel once per proline, so log2 H/L falls with PROLINE COUNT.
    The dose slope is what makes this specific -- a flat offset is incomplete labeling instead.
    Direct route: re-search the pilot with Pro6 variable and take I(Pro6)/(I(Pro6)+I(Pro0)).'''
    r = np.log2(pd.to_numeric(peptides[ratio], errors='coerce'))
    npro = peptides[seq].astype(str).str.count('P')
    ok = np.isfinite(r)
    r, npro = r[ok], npro[ok]
    slope = float(np.polyfit(npro, r, 1)[0]) if npro.nunique() > 1 else float('nan')
    return {'n_pro_free': int((npro == 0).sum()), 'n_pro_containing': int((npro > 0).sum()),
            'median_log2_HL_pro_free': round(float(r[npro == 0].median()), 4),
            'log2_HL_slope_per_proline': round(slope, 4),
            'conversion_per_proline': round(float(1 - 2 ** slope), 4)}


if __name__ == '__main__':
    import json
    import sys
    if len(sys.argv) != 2:
        sys.exit('usage: python scripts/silac_checks.py pilot_evidence.txt')
    tbl = pd.read_csv(sys.argv[1], sep='\t', low_memory=False)
    out = {'labeling': silac_labeling_efficiency(tbl)}
    if {'Sequence', 'Ratio H/L'} <= set(tbl.columns):
        out['arg_to_pro'] = arg_to_pro_shift(tbl)
    print(json.dumps(out, indent=2))
