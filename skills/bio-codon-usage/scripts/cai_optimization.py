'''Build a CAI index from highly expressed reference genes, score a query, and codon-optimize it.

CodonAdaptationIndex replaces the Bio.SeqUtils.CodonUsage module removed in Biopython 1.82.
CAI is only meaningful against an expression-biased reference (highly expressed genes of the
target host), so the index is built here from a small reference set, not a bundled table.
'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex
from Bio.Data.CodonTable import standard_dna_table

from codon_utils import calculate_cai_guarded, optimize_dna_table_aware, validate_cds

# Stand-in for highly expressed E. coli genes (ribosomal proteins, EFs), which strongly
# prefer GCT/CTG/GAC/AAA/GGT. In practice these are parsed from a FASTA with SeqIO.parse.
reference_genes = [
    Seq('ATGGCTGCTGCTCTGCTGCTGGACGACAAAAAAGGTGGTGCTCTGGACAAAGGTGCTTAA'),
    Seq('ATGGCTCTGCTGGACAAAGGTGCTGCTCTGCTGGACGACAAAAAAGGTGGTGCTCTGTAA'),
    Seq('ATGCTGGCTGACGACAAAGGTGGTGCTGCTCTGCTGAAAAAAGACGGTGCTCTGGCTTAA'),
]

def main():
    for reference in reference_genes:
        validate_cds(reference, table_id=1, policy='strict')
    cai = CodonAdaptationIndex(reference_genes, table=standard_dna_table)
    print(f'Index built from {len(reference_genes)} reference genes ({len(cai)} codon weights)')
    print(f'Relative adaptiveness of GCT (Ala): {cai["GCT"]:.3f}')

    # This query uses the host-dispreferred synonyms, so its CAI is low.
    query = Seq('ATGGCAGCATTAGATGATAAAGGATAA')
    query_score = calculate_cai_guarded(cai, query, table_id=1, policy='strict')
    print(f'\nQuery: {query_score.validation.sequence} ({len(query_score.validation.codons)} codons)')
    print(f'Genetic-code table: {query_score.validation.table_id}')
    print(f'Discarded input: {query_score.validation.discard_summary()}')
    print(f'CAI codons included/excluded: {query_score.included_codons}/{len(query_score.excluded_codons)}')
    print(f'Query CAI: {query_score.score:.3f}')

    # Optimize through a protein translated with the selected genetic-code table.
    # strict=False deliberately resolves equally preferred codons without raising.
    optimized, _, optimized_validation = optimize_dna_table_aware(
        cai, query, table_id=1, strict=False, policy='strict'
    )
    optimized_score = calculate_cai_guarded(cai, optimized, table_id=1, policy='strict')
    print(f'\nOptimized: {optimized}')
    print(f'Optimized CAI: {optimized_score.score:.3f}')

    # Protein preservation is necessary but does not screen ramp, structure, or cryptic elements.
    assert optimized.translate(table=1) == query.translate(table=1)
    print(f'Protein preserved: {optimized.translate(table=1)}')

    changes = [(i // 3, str(query)[i:i+3], str(optimized)[i:i+3])
               for i in range(0, len(query) - 2, 3) if str(query)[i:i+3] != str(optimized)[i:i+3]]
    print(f'\nCodons changed: {len(changes)}')
    for codon_index, old, new in changes:
        print(f'  codon {codon_index}: {old} -> {new}')


if __name__ == '__main__':
    main()
