"""Drop sgRNAs whose target-base editing (from CRISPResso2 output) is below a threshold.

Purpose: editing-efficiency filter to apply before hit calling in a base-editing screen.
Inputs:  a directory of CRISPResso_on_<sgrna> folders (each with
         Quantification_window_nucleotide_percentage_table.txt), 1-indexed target position within the
         quantification window, target base, threshold (default 0.5).
Output:  table sgrna_id, editing_pct (fraction 0-1), pass_filter.
Usage:   python filter_by_editing_efficiency.py results_dir --target-pos 5 --target-base C --threshold 0.5
Needs:   pandas.
"""
import pandas as pd


def filter_by_editing_efficiency(crispresso_outputs_dir, target_pos, target_base, efficiency_threshold=0.5):
    '''Drop sgRNAs that edit <efficiency_threshold of reads at target position.
    crispresso_outputs_dir: directory containing CRISPResso per-sample outputs.
    target_pos: 1-indexed position WITHIN THE QUANTIFICATION WINDOW, in column
    order -- CRISPResso2 does not emit a literal "Position" column.'''
    from pathlib import Path
    results = []
    for sample_dir in Path(crispresso_outputs_dir).glob('CRISPResso_on_*'):
        sgrna_id = sample_dir.name.replace('CRISPResso_on_', '')
        quant_file = sample_dir / 'Quantification_window_nucleotide_percentage_table.txt'
        if not quant_file.exists():
            continue
        # Real CRISPResso2 2.3.4 file layout (verified against actual output, not
        # assumed): rows are nucleotide identity (A/C/G/T/N/-, the index column);
        # columns are one per window position, header-labeled with the REFERENCE
        # base at that position (so headers repeat -- pandas suffixes duplicates
        # .1/.2/... -- select columns positionally, not by label). Values are
        # FRACTIONS in [0, 1], not 0-100, despite the filename.
        df = pd.read_csv(quant_file, sep='\t', index_col=0)
        # Schema check: fail loudly and specifically on drift instead of a bare
        # KeyError deep in a groupby/indexing call.
        if target_base not in df.index:
            raise ValueError(f"target_base={target_base!r} not in table rows {list(df.index)} ({quant_file}); "
                              "unexpected CRISPResso2 Quantification_window_nucleotide_percentage_table.txt schema")
        if not (1 <= target_pos <= df.shape[1]):
            raise ValueError(f"target_pos={target_pos} out of range for a {df.shape[1]}-position "
                              f"quantification window in {quant_file}")
        original_frac = df.loc[target_base].iloc[target_pos - 1]
        editing_pct = 1 - original_frac
        results.append({'sgrna_id': sgrna_id, 'editing_pct': editing_pct,
                         'pass_filter': editing_pct >= efficiency_threshold})
    return pd.DataFrame(results)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('crispresso_outputs_dir', help='directory containing CRISPResso_on_<sgrna> folders')
    ap.add_argument('--target-pos', type=int, required=True, help='1-indexed position within the quantification window')
    ap.add_argument('--target-base', default='C', help='reference base being edited (C for CBE, A for ABE)')
    ap.add_argument('--threshold', type=float, default=0.5)
    ap.add_argument('--out', help='write TSV here (default: print)')
    a = ap.parse_args()
    res = filter_by_editing_efficiency(a.crispresso_outputs_dir, a.target_pos, a.target_base, a.threshold)
    if a.out:
        res.to_csv(a.out, sep='\t', index=False)
    else:
        print(res.to_string(index=False))
