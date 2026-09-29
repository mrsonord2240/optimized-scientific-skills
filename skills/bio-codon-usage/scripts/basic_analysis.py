'''Basic codon usage analysis with shared CDS validation.'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq
from Bio.SeqUtils import GC123
from Bio.Data import CodonTable
from collections import Counter

from codon_utils import codon_frequencies, count_codons


def main():
    # Example coding sequence
    seq = Seq('ATGCGATCGATCGATCGATCGATCGATCGATCGTAA')
    counts, validation = count_codons(seq, table_id=1, policy='strict')
    freqs, _ = codon_frequencies(seq, table_id=1, policy='strict')
    validated_seq = Seq(validation.sequence)
    print(f'Sequence: {validated_seq}')
    print(f'Length: {len(validated_seq)} bp ({len(validation.codons)} codons)')
    print(f'Genetic-code table: {validation.table_id}')
    print(f'Discarded input: {validation.discard_summary()}')

    # Count codons
    print('\n=== Codon Counts ===')
    for codon, count in sorted(counts.items()):
        print(f'{codon}: {count}')

    # Frequencies
    print('\n=== Codon Frequencies ===')
    for codon, freq in sorted(freqs.items(), key=lambda x: x[1], reverse=True):
        bar = '#' * int(freq * 50)
        print(f'{codon}: {freq:.3f} {bar}')

    # GC at codon positions
    print('\n=== GC at Codon Positions ===')
    gc_total, gc_pos1, gc_pos2, gc_pos3 = GC123(validated_seq)
    print(f'Total GC:        {gc_total:.1f}%')
    print(f'1st position:    {gc_pos1:.1f}%')
    print(f'2nd position:    {gc_pos2:.1f}%')
    print(f'3rd position:    {gc_pos3:.1f}% (wobble)')

    # Amino acid usage
    print('\n=== Amino Acid Usage ===')
    table = CodonTable.unambiguous_dna_by_id[validation.table_id]
    aa_counts = Counter()
    for codon, count in counts.items():
        if codon in table.stop_codons:
            aa_counts['*'] += count
        elif codon in table.forward_table:
            aa_counts[table.forward_table[codon]] += count

    for aa, count in sorted(aa_counts.items(), key=lambda x: x[1], reverse=True):
        print(f'{aa}: {count}')


if __name__ == '__main__':
    main()
