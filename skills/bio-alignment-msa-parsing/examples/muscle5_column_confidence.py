'''Mask unreliable alignment columns with MUSCLE5 ensemble column confidence (CC).

Replaces GUIDANCE2, whose stand-alone package can no longer be downloaded. Pipeline (MUSCLE 5.3):

    muscle -align seqs.fa -stratified -output ens.efa          # ensemble of replicate alignments
    muscle -maxcc ens.efa -output maxcc.afa                    # stderr ends "best <name>", e.g. acb.2
    muscle -addconfseq ens.efa -output ens_cc.efa              # adds two digit-rows (CC) per replicate
    python muscle5_column_confidence.py ens_cc.efa acb.2 0.9 masked.fa

`-addconfseq` appends two sequences named `_conf_` and `_conf_2` to every replicate: the tens and
units digits of the column confidence, so "73" is CC 0.73 and "++" is 1.0. The replicate name is the
`<name` line in the .efa. Columns with CC < min_cc are removed. There is no calibrated cut-off:
run once without min_cc to see how many columns survive at each level and choose from that.
'''
# Reference: biopython 1.83+ | Verify API if version differs

import sys
from io import StringIO

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

from msa_utils import select_columns


def read_replicate(efa_path, name):
    '''(alignment, cc) for replicate `name` of a `muscle -addconfseq` .efa file; cc is one float per column.'''
    blocks, current = {}, None
    with open(efa_path) as handle:
        for line in handle:
            if line.startswith('<'):
                current = line[1:].strip()
                blocks[current] = []
            elif current is not None:
                blocks[current].append(line)
    if name not in blocks:
        raise ValueError(f'replicate {name!r} not in {efa_path}; found {sorted(blocks)}')
    records = list(SeqIO.parse(StringIO(''.join(blocks[name])), 'fasta'))
    digits = {r.id: str(r.seq) for r in records if r.id.startswith('_conf_')}
    if set(digits) != {'_conf_', '_conf_2'}:
        raise ValueError('no _conf_ rows: run `muscle -addconfseq ens.efa -output ens_cc.efa` first')
    seqs = [r for r in records if not r.id.startswith('_conf_')]

    def digit(char):
        return 10 if char == '+' else int(char)
    cc = [min(1.0, (10 * digit(t) + digit(u)) / 100) for t, u in zip(digits['_conf_'], digits['_conf_2'])]
    alignment = MultipleSeqAlignment(seqs)
    if len(cc) != alignment.get_alignment_length():
        raise ValueError('confidence rows and alignment differ in length')
    return alignment, cc


def mask_by_confidence(alignment, cc, min_cc):
    keep = [i for i, value in enumerate(cc) if value >= min_cc]
    return select_columns(alignment, keep)


if __name__ == '__main__':
    efa, name = sys.argv[1], sys.argv[2]
    alignment, cc = read_replicate(efa, name)
    print(f'{len(alignment)} sequences x {len(cc)} columns; mean CC {sum(cc) / len(cc):.3f}')
    for cutoff in (0.5, 0.7, 0.9, 0.99):
        print(f'  CC >= {cutoff}: {sum(v >= cutoff for v in cc)} columns')
    if len(sys.argv) > 4:
        masked = mask_by_confidence(alignment, cc, float(sys.argv[3]))
        SeqIO.write(masked, sys.argv[4], 'fasta')
        print(f'Wrote {masked.get_alignment_length()} columns to {sys.argv[4]}')
