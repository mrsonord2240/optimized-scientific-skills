'''Count substitutions observed in alignment (and Ti/Tv for nucleotide alignments)'''
# Reference: biopython 1.83+ (checked on 1.88) | Verify API if version differs

from collections import defaultdict

from msa_utils import example_path, is_nucleotide, load_alignment

TRANSITIONS = ({'A', 'G'}, {'C', 'T'})


def substitution_counts(alignment):
    '''Unordered mismatching residue pairs over all sequence pairs per column; input must be normalised.'''
    counts = defaultdict(int)
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        chars = [c for c in column if c != '-']
        for i, c1 in enumerate(chars):
            for c2 in chars[i+1:]:
                if c1 != c2:
                    pair = tuple(sorted([c1, c2]))
                    counts[pair] += 1
    return dict(counts)


def transition_transversion(counts):
    '''(transitions, transversions, other) from substitution_counts of a nucleotide alignment (RNA U mapped
    to T by load_alignment). Pairs involving N or other ambiguity codes are 'other'.'''
    ti = sum(v for k, v in counts.items() if set(k) in TRANSITIONS)
    tv = sum(v for k, v in counts.items() if set(k) <= set('ACGT') and set(k) not in TRANSITIONS)
    return ti, tv, sum(counts.values()) - ti - tv


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_dna.fasta'))
    print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

    subs = substitution_counts(alignment)
    total_subs = sum(subs.values())

    print(f'Total substitutions observed: {total_subs}\n')
    print('Substitution counts (sorted by frequency):')
    for pair, count in sorted(subs.items(), key=lambda x: -x[1]):
        pct = count / total_subs * 100
        print(f'  {pair[0]} <-> {pair[1]}: {count:5d} ({pct:5.1f}%)')

    if is_nucleotide(alignment):
        transitions, transversions, other = transition_transversion(subs)
        print(f'\nTransitions: {transitions}')
        print(f'Transversions: {transversions}')
        print(f'Ambiguity-code pairs (excluded from Ti/Tv): {other}')
        if transversions > 0:
            print(f'Ti/Tv ratio: {transitions/transversions:.2f}')
    else:
        print('\nProtein alignment: Ti/Tv not reported (it is defined for nucleotides only)')
