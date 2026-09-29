'''RSCU (Relative Synonymous Codon Usage) analysis'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq
from Bio.Data import CodonTable

from codon_utils import count_codons

def calculate_rscu(seq, table_id=1, policy='strict'):
    table = CodonTable.unambiguous_dna_by_id[table_id]
    counts, validation = count_codons(seq, table_id=table_id, policy=policy)

    # Build back table: amino acid -> list of codons
    back_table = {}
    for codon, aa in table.forward_table.items():
        back_table.setdefault(aa, []).append(codon)

    rscu = {}
    for aa, codons in back_table.items():
        total = sum(counts.get(c, 0) for c in codons)
        n_synonymous = len(codons)
        expected = total / n_synonymous if n_synonymous > 0 else 0
        for codon in codons:
            observed = counts.get(codon, 0)
            rscu[codon] = observed / expected if expected > 0 else 0
    return rscu, validation

def main():
    # Example sequence
    seq = Seq('ATGCTTCTGCTACTGCTTCTACTGCTGCTACTGTAA')
    rscu, validation = calculate_rscu(seq)
    print(f'Sequence: {validation.sequence}')
    print(f'Genetic-code table: {validation.table_id}')
    print(f'Discarded input: {validation.discard_summary()}')

    # Calculate RSCU
    print('\n=== RSCU Analysis ===')

    # Get codon table for amino acid names
    table = CodonTable.unambiguous_dna_by_id[validation.table_id]

    # Group by amino acid
    back_table = {}
    for codon, aa in table.forward_table.items():
        back_table.setdefault(aa, []).append(codon)

    print('RSCU values (1.0 = no bias, >1 = preferred, <1 = avoided):')
    print()

    for aa in sorted(back_table.keys()):
        codons = back_table[aa]
        if len(codons) == 1:
            continue  # Skip amino acids with only one codon

        codon_info = [(c, rscu.get(c, 0)) for c in codons]
        non_zero = [c for c, r in codon_info if r > 0]
        if not non_zero:
            continue

        print(f'{aa}:')
        for codon, r in sorted(codon_info, key=lambda x: x[1], reverse=True):
            if r > 0:
                bar = '+' * int(r * 5) if r > 1 else '-' * int((1 - r) * 5)
                status = 'preferred' if r > 1.2 else 'avoided' if r < 0.8 else 'neutral'
                print(f'  {codon}: {r:.2f} {bar} ({status})')
        print()

    # Find rare codons
    print('=== Rare Codons (RSCU < 0.5) ===')
    rare = {c: r for c, r in rscu.items() if 0 < r < 0.5}
    counts, _ = count_codons(seq, table_id=validation.table_id, policy='strict')
    for codon, r in sorted(rare.items(), key=lambda x: x[1]):
        aa = table.forward_table.get(codon, '?')
        print(f'{codon} ({aa}): RSCU={r:.2f}, count={counts.get(codon, 0)}')


if __name__ == '__main__':
    main()
