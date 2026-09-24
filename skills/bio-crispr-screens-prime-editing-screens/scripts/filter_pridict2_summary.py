# Purpose: filter a PRIDICT2 batch summary to efficient pegRNAs and keep the top N per intended edit.
# Input:   the <timestamp>_summary_<K562|HEK>_batch_summary.csv written by
#          `pridict2_pegRNA_design.py batch ... --summarize K562` (needs sequence_name and the score column).
# Output:  a CSV of the kept pegRNAs (default peg_library_filtered.csv); prints how many passed.
# Usage:   python filter_pridict2_summary.py predictions/<timestamp>_summary_K562_batch_summary.csv \
#              [--threshold 50] [--top 3] [--column PRIDICT2_0_editing_Score_deep_K562] [--out peg_library_filtered.csv]
# Checked with PRIDICT2 git HEAD 2026-09-16 output columns, pandas 3.0.
import argparse
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('summary_csv')
ap.add_argument('--threshold', type=float, default=50.0,
                help='keep pegRNAs scoring above this (0-100 scale; a project-chosen cutoff, PRIDICT2 prescribes none)')
ap.add_argument('--top', type=int, default=3, help='pegRNAs kept per sequence_name')
ap.add_argument('--column', default='PRIDICT2_0_editing_Score_deep_K562',
                help='score column; use ..._HEK for a --summarize HEK run')
ap.add_argument('--out', default='peg_library_filtered.csv')
args = ap.parse_args()

predictions = pd.read_csv(args.summary_csv)
# An empty summary (just "") means zero designs succeeded; see SKILL.md Failure Modes.
if predictions.empty or args.column not in predictions.columns:
    raise SystemExit(f'{args.summary_csv} has no {args.column} column or no rows: the batch produced no '
                     'designs (see SKILL.md Failure Modes, "Batch run exits 0 with an empty summary file")')

filtered = predictions[predictions[args.column] > args.threshold]
print(f'pegRNAs passing PRIDICT2 >{args.threshold:g}: {len(filtered)} / {len(predictions)}')

top = (filtered.sort_values(['sequence_name', args.column], ascending=[True, False])
               .groupby('sequence_name').head(args.top))
top.to_csv(args.out, index=False)
print(f'wrote {len(top)} pegRNAs to {args.out}')
