'''Shared helpers for the msa-parsing examples: gap/case normalisation, column selection, loading.

HMMER/Stockholm/A2M write '.' for gaps and soft-masked DNA is lowercase. Every helper that
compares against '-' or counts letters must see a normalised alignment, otherwise it returns
silently wrong numbers (checked on the Pfam PF00042 seed, which uses '.' gaps).
'''
# Reference: biopython 1.83+ (checked on 1.88), numpy 1.26+ | Verify API if version differs

import os
import sys

import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
STOCKHOLM_EXTENSIONS = ('.sto', '.stk', '.stockholm')


def example_path(default_name):
    '''argv[1] when given, otherwise the tiny alignment shipped in examples/data.'''
    return sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA_DIR, default_name)


def select_columns(alignment, keep, upper=None):
    '''New alignment holding only column indices `keep`.

    Copies id/name/description, per-record annotations, letter_annotations, alignment
    annotations and column_annotations (sliced to the kept columns), so a Stockholm write
    afterwards keeps its GC/GR lines. upper=True/False also maps '.' -> '-' and sets the case
    (upper-case when True, unchanged when False); upper=None leaves the sequences untouched.
    '''
    keep = list(keep)
    full = keep == list(range(alignment.get_alignment_length()))
    records = []
    for record in alignment:
        text = str(record.seq)  # once per record: str(record.seq) inside the join would be quadratic
        seq = text if full else ''.join(text[i] for i in keep)
        if upper is not None:
            seq = (seq.upper() if upper else seq).replace('.', '-')
        new = SeqRecord(Seq(seq), id=record.id, name=record.name, description=record.description,
                        dbxrefs=list(record.dbxrefs), annotations=dict(record.annotations))
        for key, values in record.letter_annotations.items():
            picked = [values[i] for i in keep]
            new.letter_annotations[key] = ''.join(picked) if isinstance(values, str) else picked
        records.append(new)
    column_annotations = {}
    for key, values in getattr(alignment, 'column_annotations', {}).items():
        picked = [values[i] for i in keep]
        column_annotations[key] = ''.join(picked) if isinstance(values, str) else picked
    return MultipleSeqAlignment(records, annotations=dict(getattr(alignment, 'annotations', {})),
                                column_annotations=column_annotations)


def normalize_alignment(alignment, upper=True):
    '''Copy of `alignment` with '.' gaps turned into '-' and (upper=True) letters upper-cased.

    Use upper=False for A2M/A3M, where lowercase marks insert states.
    '''
    return select_columns(alignment, range(alignment.get_alignment_length()), upper=upper)


def is_nucleotide(alignment, min_fraction=0.9):
    '''True when >= min_fraction of the non-gap characters are A/C/G/T/U/N (any case).'''
    text = ''.join(str(r.seq) for r in alignment).upper().replace('-', '').replace('.', '')
    if not text:
        return False
    return sum(text.count(c) for c in 'ACGTUN') / len(text) >= min_fraction


def to_array(alignment):
    '''(n_seqs, n_cols) array of characters, gaps as '-', upper case.'''
    return np.array([list(str(r.seq).upper().replace('.', '-')) for r in alignment])


def guess_format(path):
    ''''stockholm' for .sto/.stk/.stockholm files, otherwise 'fasta'.'''
    return 'stockholm' if path.lower().endswith(STOCKHOLM_EXTENSIONS) else 'fasta'


def load_alignment(path, fmt=None):
    '''Read and normalise an alignment; format is guessed from the extension when fmt is None.'''
    return normalize_alignment(AlignIO.read(path, fmt or guess_format(path)))
