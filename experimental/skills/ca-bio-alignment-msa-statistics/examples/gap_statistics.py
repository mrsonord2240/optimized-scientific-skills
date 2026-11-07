'''Analyze gap distribution in alignment'''
# Reference: biopython 1.83+ (checked on 1.88) | Verify API if version differs

from msa_utils import example_path, load_alignment


def gap_statistics(alignment):
    '''Gap counts of a normalised alignment: per sequence and per column, totals, gap-free and all-gap columns,
    and the gappiest sequence id and column index (first one on ties).'''
    num_seqs = len(alignment)
    num_cols = alignment.get_alignment_length()
    gaps_per_seq = [(r.id, str(r.seq).count('-')) for r in alignment]
    gaps_per_col = [alignment[:, i].count('-') for i in range(num_cols)]
    total_gaps = sum(g for _, g in gaps_per_seq)
    total_positions = num_seqs * num_cols
    return {
        'num_seqs': num_seqs, 'num_cols': num_cols,
        'gaps_per_seq': gaps_per_seq, 'gaps_per_col': gaps_per_col,
        'total_gaps': total_gaps, 'total_positions': total_positions,
        'gap_fraction': total_gaps / total_positions if total_positions else float('nan'),
        'gap_free_columns': sum(1 for g in gaps_per_col if g == 0),
        'all_gap_columns': sum(1 for g in gaps_per_col if g == num_seqs),
        'gappiest_sequence': max(gaps_per_seq, key=lambda x: x[1])[0] if gaps_per_seq else None,
        'gappiest_column': max(range(num_cols), key=lambda i: gaps_per_col[i]) if num_cols else None,
    }


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_protein.fasta'))
    stats = gap_statistics(alignment)
    num_seqs, num_cols = stats['num_seqs'], stats['num_cols']
    gaps_per_seq, gaps_per_col = stats['gaps_per_seq'], stats['gaps_per_col']
    total_gaps, total_positions = stats['total_gaps'], stats['total_positions']

    print(f'Alignment: {num_seqs} sequences, {num_cols} columns\n')

    print('Gap Statistics:')
    print(f'  Total gaps: {total_gaps}')
    print(f'  Total positions: {total_positions}')
    print(f'  Overall gap fraction: {total_gaps/total_positions*100:.1f}%')

    print('\nGaps per sequence:')
    for seq_id, gaps in sorted(gaps_per_seq, key=lambda x: -x[1]):
        pct = gaps / num_cols * 100
        print(f'  {seq_id}: {gaps} ({pct:.1f}%)')

    gap_free, all_gap = stats['gap_free_columns'], stats['all_gap_columns']
    print(f'\nColumn statistics:')
    print(f'  Gap-free columns: {gap_free} ({gap_free/num_cols*100:.1f}%)')
    print(f'  All-gap columns: {all_gap}')

    if any(g > 0 for g in gaps_per_col):
        # 0.5 (50%) threshold: columns with majority gaps often indicate alignment uncertainty
        # or indel events. Lower threshold (0.3) for stringent analysis; higher (0.7) if gaps expected.
        gappy_cols = [(i, g) for i, g in enumerate(gaps_per_col) if g > num_seqs * 0.5]
        if gappy_cols:
            print(f'\nColumns with >50% gaps: {len(gappy_cols)}')
            for i, g in gappy_cols[:10]:
                print(f'  Column {i}: {g}/{num_seqs} gaps')
