# Reference: PRIDICT2 (uzh-dqbm-cmi/PRIDICT2, git HEAD 2026-09-16), pandas 2.2+ | Verify API if version differs
# PBS/RTT extraction uses a hand-rolled revcomp (no biopython dependency needed for this script).
#
# pegRNA library design with PRIDICT2 efficiency prediction.
# PRIDICT2 is invoked via CLI (pridict2_pegRNA_design.py batch).
# Filters to high-efficiency pegRNAs before library synthesis.

import pandas as pd
import subprocess
from pathlib import Path
import re

# === INPUTS ===
# Intended variants: variant_id, chrom, pos, ref, alt, context
variants_df = pd.read_csv('intended_variants.csv')

# === PRIDICT2 CLI invocation ===
# Build the batch input CSV expected by PRIDICT2: sequence_name, editseq (with (REF/ALT) notation)
# Then call (CSV under ./input/, out/ must exist, and no .csv already in out/):
#   mkdir -p input/ out/; mv variants.csv input/
#   python pridict2_pegRNA_design.py batch --input-fname variants.csv --output-dir out/ --cores 3 --summarize K562
# This example uses a placeholder predict_pridict2() that simulates the CLI output.

# === STEP 1: GENERATE pegRNA CANDIDATES ===
# PBS/RTT geometry (verified against real CRISPResso2 2.3.4 output and the real
# PRIDICT2 CLI's own PBSrevcomp/RTrevcomp columns -- see SKILL.md's pegRNA
# Architecture section):
#   - The nick is 3 nt upstream of the PAM (cut_pos = pam_pos - 3).
#   - PBS is the region UPSTREAM of the nick, i.e. part of the protospacer
#     itself -- it must never include PAM bases.
#   - RTT is the region DOWNSTREAM of the nick (through the PAM and beyond) --
#     this is the window PE actually replaces, and where the edit lives
#     (1-30 nt from the nick).
#   - Both PBS and RTT are reverse-complemented relative to the protospacer
#     strand to get the pegRNA-sense sequence, and the 3' extension is
#     RTT-then-PBS (not PBS-then-RTT).
_COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def _revcomp(s):
    return ''.join(_COMP[b] for b in reversed(s))


# Try the PRIDICT2-typical lengths first (13 nt PBS, 18 nt RTT), then widen.
PBS_LENGTHS = (13, 12, 11, 14, 15, 10, 9, 8)
RTT_LENGTHS = (18, 15, 20, 13, 25, 10, 30)


def find_pegrna_candidates(chrom, pos, ref, alt, context_seq, edit_position_in_context=30):
    '''Find pegRNAs that can install this variant.
    `context_seq` must place the edit at `edit_position_in_context` (default 30,
    matching this script's 60-nt context convention -- see variants_to_design.csv).
    Returns list of dicts: spacer, pbs, rtt, extension_seq (5'->3', RTT then PBS),
    pbs_length, rtt_length, pam_strand, edit_dist. Empty list if no NGG PAM gives
    a workable PBS/RTT window -- callers must check for this (see Step 3).'''
    candidates = []
    for strand in ['+', '-']:
        if strand == '+':
            seq, edit_pos, alt_on_strand = context_seq, edit_position_in_context, alt
        else:
            # Search the other strand in its own 5'->3' coordinate frame; the
            # edit position and installed base both have to flip with it.
            seq = _revcomp(context_seq)
            edit_pos = len(context_seq) - 1 - edit_position_in_context
            alt_on_strand = _COMP[alt]
        for pam_match in re.finditer(r'(?=([ACGT]GG))', seq):
            pam_pos = pam_match.start()
            if pam_pos < 20:
                continue  # no room for a full 20nt spacer upstream of this PAM
            cut_pos = pam_pos - 3
            edit_dist = edit_pos - cut_pos  # nt from nick to edit (must be downstream)
            if not (1 <= edit_dist <= 30):
                continue
            spacer = seq[pam_pos - 20:pam_pos]
            # Find the first (pbs_length, rtt_length) pair that fits this locus.
            found = None
            for pbs_length in PBS_LENGTHS:
                if cut_pos - pbs_length < 0:
                    continue
                for rtt_length in RTT_LENGTHS:
                    if edit_dist >= rtt_length or cut_pos + rtt_length > len(seq):
                        continue
                    found = (pbs_length, rtt_length)
                    break
                if found:
                    break
            if not found:
                continue
            pbs_length, rtt_length = found
            pbs_geno = seq[cut_pos - pbs_length:cut_pos]
            rtt_geno = list(seq[cut_pos:cut_pos + rtt_length])
            rtt_geno[edit_dist] = alt_on_strand
            pbs = _revcomp(pbs_geno)
            rtt = _revcomp(''.join(rtt_geno))
            candidates.append({
                'spacer': spacer,
                'pbs': pbs,
                'rtt': rtt,
                'pbs_length': pbs_length,
                'rtt_length': rtt_length,
                'pam_strand': strand,
                'edit_dist': edit_dist,
                'extension_seq': rtt + pbs,  # RTT then PBS, not PBS then RTT
            })
    return candidates

# === STEP 2: PRIDICT2 PREDICTION ===
# In production: write candidates to a CSV with PRIDICT2's batch format
# (sequence_name, editseq with (REF/ALT) notation) and call:
#   mkdir -p input/ predictions/; mv candidates.csv input/
#   python pridict2_pegRNA_design.py batch --input-fname candidates.csv \
#       --output-dir predictions/ --cores 3 --summarize K562
# Note: --output-dir must be empty of .csv files before running with --summarize.
# Then read predictions/<timestamp>_summary_K562_batch_summary.csv. The placeholder below simulates.
def predict_pridict2(spacer, pbs, rtt, context):
    '''Placeholder for PRIDICT2 batch CLI output parsing.
    In real use: parse predictions/<timestamp>_summary_K562_batch_summary.csv after running PRIDICT2 CLI.'''
    pbs_gc = (pbs.count('G') + pbs.count('C')) / len(pbs)
    sim_efficiency = max(0.05, min(0.95, 0.4 + 0.3 * (0.5 - abs(pbs_gc - 0.45))))
    return {'efficiency': sim_efficiency, 'indel_rate': 0.02, 'scaffold_incorp': 0.03}

# === STEP 3: SCREEN ALL CANDIDATES ===
all_pegrnas = []
variants_with_no_pam = []
for _, var in variants_df.iterrows():
    candidates = find_pegrna_candidates(var['chrom'], var['pos'],
                                          var['ref'], var['alt'], var['context'])
    if not candidates:
        # No NGG PAM gave a workable PBS/RTT window for this variant -- not a
        # bug, a real PE-design dead end (see Failure Modes: "Library missing
        # intended variant"). Skip it, don't crash the whole batch.
        variants_with_no_pam.append(var['variant_id'])
        continue
    for c in candidates:
        pred = predict_pridict2(c['spacer'], c['pbs'], c['rtt'], var['context'])
        all_pegrnas.append({
            'variant_id': var['variant_id'],
            'spacer': c['spacer'],
            'pbs': c['pbs'],
            'rtt': c['rtt'],
            'extension_seq': c['extension_seq'],
            'pam_strand': c['pam_strand'],
            'edit_dist': c['edit_dist'],
            'pridict2_efficiency': pred['efficiency'],
            'pridict2_indel': pred['indel_rate'],
            'pridict2_scaffold': pred['scaffold_incorp'],
        })

if variants_with_no_pam:
    print(f'No PE-designable PAM/PBS/RTT window for {len(variants_with_no_pam)} variant(s): '
          f'{variants_with_no_pam} -- excluded before PRIDICT2 scoring, not a script error.')

if not all_pegrnas:
    raise ValueError(
        'No pegRNA candidates found for any variant in this batch -- check that '
        'each context places the edit at edit_position_in_context (default 30) '
        'and has an NGG PAM within 30 nt of it. This is a design dead end, not '
        'a code bug; see SKILL.md Failure Modes: "Library missing intended variant".'
    )

pegrnas_df = pd.DataFrame(all_pegrnas)

# === STEP 4: FILTER AND RANK ===
# Project-chosen filter: PRIDICT2 efficiency >50% (the paper prescribes no numeric cutoff)
efficient = pegrnas_df[pegrnas_df['pridict2_efficiency'] > 0.50]

# Pick top 3 per variant for library
top3 = (efficient.sort_values(['variant_id', 'pridict2_efficiency'],
                                ascending=[True, False])
                  .groupby('variant_id')
                  .head(3))

# === STEP 5: REPORT ===
n_variants_total = variants_df['variant_id'].nunique()
n_variants_covered = top3['variant_id'].nunique()
print(f'Variants attempted: {n_variants_total}')
print(f'Variants with >=1 pegRNA passing PRIDICT2 >50%: {n_variants_covered}')
print(f'Total pegRNAs (top 3 per variant): {len(top3)}')

# Flag variants without coverage
uncovered = set(variants_df['variant_id']) - set(top3['variant_id'])
print(f'Uncovered variants requiring alternative (BE, SpRY-PE, or accept exclusion): {len(uncovered)}')

# === STEP 6: EXPORT ===
top3.to_csv('peg_library_filtered.csv', index=False)
if uncovered:
    pd.DataFrame({'variant_id': list(uncovered)}).to_csv('uncovered_variants.csv', index=False)

print('\nNext steps:')
print('  1. Pilot top pegRNAs at chromatin-accessible loci first')
print('  2. Synthesize library via oligo pool')
print('  3. Run pooled PE screen at MOI 0.3 in PE-validated cell line')
print('  4. Endpoint amplicon sequencing + CRISPResso2 PE mode')
print('  5. Filter to >5% intended-edit pegRNAs')
print('  6. MAGeCK MLE per-variant hit calling')
