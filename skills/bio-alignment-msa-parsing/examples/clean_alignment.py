'''Clean alignment by removing gappy columns and sequences'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from msa_utils import example_path, guess_format, normalize_alignment, select_columns


def find_gappy_columns(alignment, threshold=0.5):
    alignment = normalize_alignment(alignment)
    num_seqs = len(alignment)
    return [i for i in range(alignment.get_alignment_length())
            if alignment[:, i].count('-') / num_seqs >= threshold]


def remove_gappy_columns(alignment, threshold=0.5):
    gappy = set(find_gappy_columns(alignment, threshold))
    keep = [i for i in range(alignment.get_alignment_length()) if i not in gappy]
    return select_columns(alignment, keep)  # keeps record and column annotations


def filter_by_gap_content(alignment, max_gap_fraction=0.2):
    normalized = normalize_alignment(alignment)
    fractions = [str(r.seq).count('-') / len(r.seq) for r in normalized]
    kept = [r for r, f in zip(alignment, fractions) if f <= max_gap_fraction]
    if not kept:
        raise ValueError(f'max_gap_fraction={max_gap_fraction} removes all {len(alignment)} sequences '
                         f'(lowest gap fraction is {min(fractions):.2f})')
    return MultipleSeqAlignment(kept, annotations=alignment.annotations,
                                column_annotations=alignment.column_annotations)


if __name__ == '__main__':
    path = example_path('example_alignment.fasta')
    alignment = AlignIO.read(path, guess_format(path))
    print(f'Original: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

    # threshold=0.5: Remove columns with >=50% gaps. Standard cutoff for phylogenetics.
    # Use 0.3 for stringent filtering; 0.7 if expecting many indels.
    cleaned = remove_gappy_columns(alignment, threshold=0.5)
    print(f'After column cleaning: {len(cleaned)} sequences, {cleaned.get_alignment_length()} columns')

    # max_gap_fraction=0.2: Remove sequences with >20% gaps. Excludes fragmentary sequences.
    # Use 0.1 for high-quality datasets; 0.3-0.4 for diverse or ancient sequences.
    cleaned = filter_by_gap_content(cleaned, max_gap_fraction=0.2)
    print(f'After sequence filtering: {len(cleaned)} sequences, {cleaned.get_alignment_length()} columns')

    AlignIO.write(normalize_alignment(cleaned, upper=False), 'cleaned_alignment.fasta', 'fasta')  # FASTA gaps as '-'
    print('Saved to cleaned_alignment.fasta')
