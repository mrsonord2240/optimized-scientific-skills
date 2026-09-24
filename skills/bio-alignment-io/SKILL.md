---
name: bio-alignment-io
category: Data Analysis
description: Read, write, and convert multiple sequence alignment files using Biopython Bio.AlignIO. Supports Clustal, PHYLIP, Stockholm, FASTA, Nexus, and other alignment formats for phylogenetics and conservation analysis. Use when reading, writing, or converting alignment file formats.
tool_type: python
primary_tool: Bio.AlignIO
license: MIT
author: GPTomics
---

## Version Compatibility

Checked on Biopython 1.88 and pyhmmer 0.12.3 (2026-09-19), downstream tools on RAxML-NG 1.2.2 / 2.0.3, PhyML 3.3.20220408 / 3.3.20260528, MrBayes 3.2.7, IQ-TREE 3.1.3 (2026-09-20); patterns need Biopython 1.83+.

```bash
pip install biopython
pip install pyhmmer   # optional: Pfam-scale streaming, A2M reading
```

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment File I/O

Read, write, and convert multiple sequence alignment files in various formats.

## Required Import

**Goal:** Load modules for reading, writing, and manipulating multiple sequence alignments.

**Approach:** Import AlignIO for file I/O and supporting classes for programmatic alignment construction.

```python
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
```

## Format Coverage Map

Three Python libraries cover the alignment-format space, with overlapping but non-identical support. Pick by what is actually required. Checked on Biopython 1.88 / pyhmmer 0.12.3.

| Format | `Bio.AlignIO` (`MultipleSeqAlignment`) | `Bio.Align` (`Alignment`) | `pyhmmer.easel` (`format=`) | Notes |
|--------|---------------|----------------------|-----------------|-------|
| Aligned FASTA | R/W | R/W | R/W (`afa`) | Most portable; loses annotations |
| Clustal | R/W | R/W | R/W | Clustal conservation marks NOT round-tripped |
| PHYLIP (interleaved/sequential/relaxed) | R/W | R/W (strict 10-char names only: no relaxed variant; long ids raise on read and are silently truncated on write) | R/W (`phylip`, `phylips`) | See PHYLIP pitfalls |
| Stockholm | R/W | R (plain files only; TypeError on real Pfam), W (AttributeError) | R/W | Only format preserving GS/GR/GC annotations (GF header lines are dropped by AlignIO) |
| NEXUS | R/W | R/W | -- | MrBayes / PAUP* input; write needs a molecule type |
| MAF (Multiple Alignment Format) | R/W | R/W | -- | UCSC whole-genome alignments |
| A2M | -- (`'fasta'` parser works on padded A2M only) | R (padded A2M), W (needs `column_annotations['state']`, so only for alignments read from A2M) | R/W | HMMER `hmmalign` writes ragged A2M |
| A3M | -- | -- | -- | HH-suite / ColabFold; convert to A2M first (below) |
| MSF (GCG) | R | R | -- | GCG legacy |
| EMBOSS / Mauve XMFA / FASTA-m10 | R (Mauve also W; FASTA-m10 AlignIO only) | R (Mauve also W; no FASTA-m10) | -- | Mostly one-way |
| PSL / chain / BED / SAM / exonerate / bigMaf / bigPsl / bigBed | -- | R/W (pairwise alignments, not MSAs) | -- | Use Kent tools for manipulation |
| HHR / tabular (BLAST) | -- | R only (search hits, not MSAs) | -- | |

**Formats NOT in BioPython** (use dedicated tools):

| Format | Tool | Why |
|--------|------|-----|
| HAL | progressiveCactus, halTools | HDF5-backed multi-genome alignments at TB scale |
| net | UCSC Kent tools (`chainNet`) | Pairwise genome alignment |
| AXT | BLASTZ / lastz native | Pairwise alignment blocks |
| GFA / rGFA | `vg`, `odgi`, `pggb`, gfatools | Pangenome graph |
| GAF | `vg surject`, `vg call` | Graph alignment format (read-to-graph) |

Use `Bio.AlignIO` for MSA files, especially Stockholm and NEXUS; use `Bio.Align` (`Alignment` objects with `.counts()` and `.substitutions`, see "Alternative: Bio.Align Module I/O") when those features are needed. For multi-gigabyte Stockholm databases such as Pfam-A.full, read [Streaming large Stockholm databases](references/stockholm-streaming.md) before choosing an in-memory parser.

## Reference Files

Read these only for their named format or scale concern:

| Reference | Read when |
|-----------|-----------|
| [MAF coordinates](references/maf-coordinates.md) | Mapping a MAF row to plus-strand genome coordinates. |
| [A2M and A3M conventions](references/a2m-a3m.md) | Reading HMMER, HH-suite, or ColabFold alignment output. |
| [Streaming Stockholm databases](references/stockholm-streaming.md) | Iterating Pfam-A.full, BFD, or another alignment collection too large for memory. |

## Reading Alignments

**"Read an alignment file"** -> Parse an alignment file into an alignment object with sequences and metadata accessible.

**Goal:** Load alignment data from files in various formats (Clustal, PHYLIP, Stockholm, FASTA).

**Approach:** Use `AlignIO.read()` for single-alignment files or `AlignIO.parse()` for files containing multiple alignments.

### Single Alignment File
```python
from Bio import AlignIO

alignment = AlignIO.read('alignment.aln', 'clustal')
print(f'Alignment length: {alignment.get_alignment_length()}')
print(f'Number of sequences: {len(alignment)}')
```

### Multiple Alignments in One File
```python
for alignment in AlignIO.parse('multi_alignment.sto', 'stockholm'):
    print(f'Alignment with {len(alignment)} sequences, length {alignment.get_alignment_length()}')
```

### Read as List
```python
alignments = list(AlignIO.parse('alignments.phy', 'phylip'))
print(f'Read {len(alignments)} alignments')
```

## Writing Alignments

**Goal:** Save alignment data to files in standard formats for downstream tools or archival.

**Approach:** Use `AlignIO.write()` with the target format specifier, supporting single or multiple alignments and file handles.

### Write Single Alignment
```python
AlignIO.write(alignment, 'output.fasta', 'fasta')
```

### Write Multiple Alignments
```python
alignments = [alignment1, alignment2, alignment3]
count = AlignIO.write(alignments, 'output.sto', 'stockholm')
print(f'Wrote {count} alignments')
```

### Write to Handle
```python
with open('output.aln', 'w') as handle:
    AlignIO.write(alignment, handle, 'clustal')
```

## Format Conversion

**"Convert alignment format"** -> Transform an alignment file from one format to another (e.g., Clustal to PHYLIP).

**Goal:** Convert alignment files between formats for compatibility with different analysis tools.

**Approach:** Use `AlignIO.convert()` for direct one-step conversion, or read-modify-write for cases requiring intermediate manipulation.

### Direct Conversion (Most Efficient)
```python
AlignIO.convert('input.aln', 'clustal', 'output.phy', 'phylip-relaxed')
```

### NEXUS Output Needs a Molecule Type
Readers for Clustal, Stockholm, PHYLIP and FASTA leave `molecule_type` unset, and the NEXUS writer raises `ValueError: Need the molecule type to be defined` (leaving a 0-byte file). Pass it to `convert()`, or set it on each record before `write()` (`'DNA'`, `'RNA'` or `'protein'`):
```python
AlignIO.convert('input.sto', 'stockholm', 'output.nex', 'nexus', molecule_type='DNA')

alignment = AlignIO.read('input.aln', 'clustal')
for record in alignment:
    record.annotations['molecule_type'] = 'DNA'
AlignIO.write(alignment, 'output.nex', 'nexus')
```
The value is copied into `datatype=` unchecked: label a protein alignment `'DNA'` and the file says `datatype=dna` with no error, so set it from the data. `examples/convert_formats.py` recognizes the complete IUPAC nucleotide alphabet (including ambiguity codes), refuses mixed `T`/`U` data and any DNA/RNA/protein override that contradicts the residues before it writes files, then re-reads its NEXUS output.

**MrBayes needs quote-free ids.** Biopython single-quotes ids containing punctuation such as `-`, `+`, `:` or a space (Pfam/UniProt `GLB2_LUMTE/31-141`; `/` and `.` alone are not quoted), and MrBayes 3.2.7 stops with ``Instead found ''' in command 'Matrix'``. Replace every character outside `[A-Za-z0-9_]` before writing, and keep the ids unique (checked: MrBayes 3.2.7 loads the sanitized file and runs `mcmc`):
```python
import re
for record in alignment:
    record.id = re.sub(r'[^A-Za-z0-9_]', '_', record.id)
    record.name, record.description = record.id, ''
assert len({r.id for r in alignment}) == len(alignment), 'ids collide after sanitizing'
```

### Manual Conversion (When Modification Needed)
```python
alignment = AlignIO.read('input.aln', 'clustal')
# ... modify alignment ...
AlignIO.write(alignment, 'output.fasta', 'fasta')
```

## Accessing Alignment Data

**Goal:** Navigate and extract data from alignment objects including sequences, columns, and slices.

**Approach:** Use iteration, indexing, and column slicing on the alignment object.

```python
alignment = AlignIO.read('alignment.aln', 'clustal')

# Iterate over sequences
for record in alignment:
    print(f'{record.id}: {record.seq}')

# Access by index
first_seq = alignment[0]
last_seq = alignment[-1]

# Slice columns
column_slice = alignment[:, 10:20]  # Columns 10-19

# Get specific column
column = alignment[:, 5]  # Column 5 as string

# Subset of sequences, and sequences plus columns together
subset = alignment[0:5]              # First 5 sequences
region = alignment[0:5, 50:150]      # 5 sequences, columns 50-149

seq_ids = [record.id for record in alignment]
```

Column slices past the alignment length silently return 0 columns (and 0-length records on write); check `alignment.get_alignment_length()` first.

## Creating Alignments Programmatically

**Goal:** Build an alignment object from sequences defined in code rather than read from a file.

**Approach:** Construct SeqRecord objects with gap characters and wrap them in a MultipleSeqAlignment.

```python
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

records = [
    SeqRecord(Seq('ACTGACTGACTG'), id='seq1'),
    SeqRecord(Seq('ACTGACT-ACTG'), id='seq2'),
    SeqRecord(Seq('ACTG-CTGACTG'), id='seq3'),
]
alignment = MultipleSeqAlignment(records)
AlignIO.write(alignment, 'new_alignment.fasta', 'fasta')
```

## Format Selection for Downstream Tools

Choosing the output format depends on which downstream tool consumes the alignment:

| Downstream Tool | Required Format | BioPython Format String |
|----------------|-----------------|------------------------|
| RAxML-NG, IQ-TREE | PHYLIP (relaxed) | `'phylip-relaxed'` |
| MrBayes | NEXUS (quote-free ids) | `'nexus'` |
| PAUP* | NEXUS or PHYLIP | `'nexus'` or `'phylip'` |
| HMMER, Infernal | Stockholm | `'stockholm'` |
| Pfam/Rfam databases | Stockholm | `'stockholm'` |
| PAML/codeml | PHYLIP (sequential) | `'phylip-sequential'` |
| Most tools | FASTA | `'fasta'` |

### Annotation Preservation

Not all formats support annotations. Converting between formats can silently discard metadata:

| Format | Sequence Annotations | Column Annotations | Secondary Structure |
|--------|---------------------|-------------------|-------------------|
| Stockholm | Yes (GS/GR lines) | Yes (GC lines) | Yes (SS_cons) |
| NEXUS | Partial (SETS block) | Via CHARSET | No |
| Clustal | No (conservation marks not parsed) | No | No |
| PHYLIP | No | No | No |
| FASTA | No | No | No |

Converting Stockholm to FASTA or PHYLIP discards all annotations, secondary structure markup, and per-residue quality scores. Even Stockholm to Stockholm through `AlignIO` drops the `#=GF` header lines (`ID`, `AC`, `DE`, ...; verified on Pfam PF00042). If annotations matter, keep the original Stockholm file as the master copy.

## Format-Specific Notes

### PHYLIP Format Pitfalls

PHYLIP has two incompatible variants (interleaved vs sequential) and two name-length modes (strict vs relaxed). Mixing them up usually fails loudly (`ValueError` on read), but strict names are truncated to 10 characters.

**Strict PHYLIP** truncates sequence names to exactly 10 characters. Writing distinct names that share a 10-character prefix (e.g., `Homo_sapiens_chr1` and `Homo_sapiens_chr2`) raises `ValueError: Repeated name 'Homo_sapie'` (`phylip-sequential` truncates and raises the same way). The silent case is *reading* a foreign strict file whose names collide: both records are accepted with the same id. Names that stay unique are truncated without warning.

```python
# Strict PHYLIP (10-char names, interleaved) -- only for tools requiring it
alignment = AlignIO.read('file.phy', 'phylip')

# Sequential PHYLIP (10-char names, one sequence at a time) -- PAML/codeml.
# Shorten ids to <=10 unique characters first (NCBI headers such as 'lcl|NM_001...' collide)
alignment = AlignIO.read('file.phy', 'phylip-sequential')

# Relaxed PHYLIP (no name limit) -- RAxML-NG, IQ-TREE (recommended default)
alignment = AlignIO.read('file.phy', 'phylip-relaxed')

# Always prefer phylip-relaxed for writing unless the downstream tool
# specifically requires strict format
AlignIO.write(alignment, 'output.phy', 'phylip-relaxed')
```

**Clustal IDs:** Biopython's Clustal writer truncates identifiers to 30 characters without a collision check. After writing Clustal, re-read it and assert the ids remain unique before using it as an input to another tool.

#### PHYLIP-Relaxed Dialect Mismatches Between Tree Tools

Biopython's `'phylip-relaxed'` writes a single space between name and sequence. RAxML-NG and IQ-TREE accept this (checked with names of 138 characters); PhyML and RAxML-NG reject punctuation in names (table below); PAML's codeml expects sequential format with name-truncation behaviour distinct from interleaved. Common silent failures:

| Symptom | Cause | Fix |
|---------|-------|-----|
| RAxML-NG: `ERROR: Invalid character in sequence N at position P: *` | `*` in a nucleotide alignment (RAxML-NG 1.2.2 and 2.0.3 accept `*` in protein data as undetermined, so a stop codon is not scored) | Replace with `N` or `-` |
| IQ-TREE 3.1.3: `WARNING: Some sequence names are changed` | Foreign file with `:` in a name (IQ-TREE renames it; Biopython's relaxed writer already turns `:` into `|` and drops `(` `,`) | Sanitize names yourself: `re.sub(r'[^A-Za-z0-9_]', '_', record.id)` |
| PhyML: `Character ':' is not permitted in sequence name` (same for `,`) | Punctuation in a name (PhyML accepts `(` `)` but writes them into the Newick tree; 138-character names are kept intact, checked on 3.3.20220408 and 3.3.20260528) | Sanitize names with the same `[A-Za-z0-9_]`-only recipe |
| RAxML-NG: `ERROR: Following taxon name contains invalid characters` | `:` `,` `(` or `[` in a name | Sanitize names with the same `[A-Za-z0-9_]`-only recipe |
| codeml: `Error in sequence data file ... separate the sequence from its name by 2 or more spaces` | Used `phylip-relaxed` (single space, interleaved) instead of `phylip-sequential` | codeml requires sequential with short unique names |

IQ-TREE has no validate-only flag (`--check` is invalid). Before a long run, check the input with a zero-iteration run (checked on IQ-TREE 3.1.3; a wrong sequence length gives `ERROR: Line N: Sequence X has wrong sequence length`), or re-read the file with `AlignIO.read`:

```bash
iqtree3 -s file.phy -n 0 -m LG -redo -pre check   # prints 'Alignment has N sequences with M columns'
```

### MAF Block Coordinate Conventions

Read [MAF coordinates](references/maf-coordinates.md) before mapping MAF rows to a genome: minus-strand `start` values require a coordinate conversion and Biopython reports `strand` as an integer.

### Stockholm Format Annotations

Stockholm format (used by Pfam, Rfam, HMMER) supports four annotation line types:

| Line Prefix | Scope | Description | Example |
|-------------|-------|-------------|---------|
| `#=GF` | File | Alignment-level metadata (ID, accession, description) | `#=GF AC PF00001` |
| `#=GC` | Column | Per-column annotation (1 char per alignment column) | `#=GC SS_cons ..(((...)))..` |
| `#=GS` | Sequence | Per-sequence free text (organism, description) | `#=GS seq1 OS Homo sapiens` |
| `#=GR` | Residue | Per-residue annotation (1 char per residue) | `#=GR seq1 SS ..HHH..EEE..` |

Common GC annotations: `SS_cons` (consensus secondary structure), `RF` (reference coordinates), `seq_cons` (consensus sequence).

**WUSS notation** in RNA `#=GC SS_cons` lines uses nested bracket pairs (`<>`, `()`, `[]`, `{}`) for paired bases and characters like `_`, `-`, `,`, `:`, `.`, `~` for unpaired regions; pseudoknots use upper/lower-case letter pairs (`Aa`, `Bb`). Consult the Infernal user guide for the full character table before writing or parsing custom SS_cons strings.

```python
alignment = AlignIO.read('pfam.sto', 'stockholm')

for record in alignment:
    print(record.id, record.annotations)
    if 'secondary_structure' in record.letter_annotations:
        print(f'  SS: {record.letter_annotations["secondary_structure"]}')

ss_cons = alignment.column_annotations.get('secondary_structure')
```

**Round-trip caveat:** re-reading and re-writing as Stockholm preserves GS/GR/GC, but dropped/added sequences invalidate the per-residue annotations -- regenerate annotations after edits.

**Pfam-style `name/start-end` identifier convention:** Pfam, Rfam, and Dfam Stockholm IDs (e.g. `Q9Y6Y0/45-198`) encode a 1-based inclusive region. `AlignIO` already copies the region into `record.annotations['start']` / `['end']` (and the name into `['accession']`) but leaves the suffix on `record.id`. IQ-TREE 3.1.3 accepts ids with `/`; strip the suffix from `record.id` only for a tool that rejects it, then restore it after.

### A2M / A3M Conventions

Read [A2M and A3M conventions](references/a2m-a3m.md) before treating HMMER, HH-suite, or ColabFold output as a rectangular alignment. It covers ragged `hmmalign` output and the first-record rule for `reformat.pl`.

### Streaming Large Stockholm Databases

For Pfam-A.full, BFD, or another collection too large for memory, use [Streaming Stockholm databases](references/stockholm-streaming.md) instead of an in-memory `AlignIO` recipe.

## Batch Processing Multiple Files

**Goal:** Convert a directory of alignment files from one format to another in bulk.

**Approach:** Glob for input files and iterate, reading each alignment and writing to the target format.

```python
from pathlib import Path

input_dir = Path('alignments/')
output_dir = Path('converted/')
output_dir.mkdir(exist_ok=True)

for input_file in input_dir.glob('*.aln'):
    alignment = AlignIO.read(input_file, 'clustal')
    output_file = output_dir / f'{input_file.stem}.fasta'
    AlignIO.write(alignment, output_file, 'fasta')
```

## Alternative: Bio.Align Module I/O

**Goal:** Use the newer Bio.Align module for alignment I/O with access to features like counts and substitutions.

**Approach:** Use `Align.read()`, `Align.parse()`, and `Align.write()` which return `Alignment` objects instead of `MultipleSeqAlignment`. Stick to FASTA, Clustal, PHYLIP and MAF here; on Biopython 1.88 Stockholm through `Bio.Align` fails (see Common Errors) and its NEXUS writer needs `molecule_type`, so use `Bio.AlignIO` for those.

```python
from Bio import Align

# Read single alignment (returns Alignment object)
alignment = Align.read('alignment.aln', 'clustal')
print(alignment.shape, alignment.counts())   # (n_seqs, n_columns), AlignmentCounts

# Parse multiple alignments (e.g. MAF blocks)
for block in Align.parse('blocks.maf', 'maf'):
    print(f'Alignment with {len(block)} sequences')

# Write alignment
Align.write(alignment, 'output.fasta', 'fasta')
```

### When to Use Which

| Use Case | Module |
|----------|--------|
| Legacy code, MultipleSeqAlignment needed | `Bio.AlignIO` |
| Stockholm (annotations), NEXUS, `molecule_type` handling | `Bio.AlignIO` |
| Modern features (counts, substitutions) on FASTA / Clustal / PHYLIP / MAF | `Bio.Align` |
| Format conversion | `Bio.AlignIO` (works for every format above) |
| Working with pairwise alignments | `Bio.Align` |

## Quick Reference: Common Operations

| Task | Code |
|------|------|
| Read single alignment | `AlignIO.read(file, format)` |
| Read multiple alignments | `AlignIO.parse(file, format)` |
| Write alignment(s) | `AlignIO.write(align, file, format)` |
| Convert format | `AlignIO.convert(in_file, in_fmt, out_file, out_fmt)` |
| Get length | `alignment.get_alignment_length()` |
| Get sequence count | `len(alignment)` |
| Slice columns | `alignment[:, start:end]` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `ValueError: No records found in handle` | Empty file | Check file path and format |
| `ValueError: More than one record found in handle` | Multiple alignments with `read()` | Use `parse()` instead |
| `ValueError: Sequences must all be the same length` | Unaligned or ragged input (e.g. `hmmalign` A2M, raw A3M) | Align first, or pad (see A2M / A3M Conventions) |
| `ValueError: Unknown format 'a2m'` | Format string not supported by `AlignIO` | Use `'fasta'` for padded A2M, or `Align.read(..., 'a2m')` |
| `ValueError: Need the molecule type to be defined` | NEXUS write without `molecule_type` | See "NEXUS Output Needs a Molecule Type" |
| `TypeError: Any per-letter annotation should be a Python sequence ...` from `Align.read/parse(..., 'stockholm')` | Bio.Align Stockholm reader on real Pfam (Biopython 1.88) | Use `AlignIO.read(..., 'stockholm')` |
| `AttributeError: ... no attribute 'column_annotations'` from `Align.write(..., 'stockholm')` | Bio.Align Stockholm writer needs an annotated alignment | Use `AlignIO.write(..., 'stockholm')` |

## Related Skills

- alignment/multiple-alignment - Run MSA tools (MAFFT, MUSCLE5, ClustalOmega) to generate alignments
- alignment/pairwise-alignment - Create pairwise alignments with PairwiseAligner
- alignment/msa-parsing - Analyze alignment content and annotations
- alignment/msa-statistics - Calculate conservation and identity
- alignment/structural-alignment - Foldseek/TM-align outputs and Foldmason `result_aa.fa` / `result_3di.fa` MSAs (FASTA-loadable; the per-column LDDT report is HTML, not BioPython-parseable)
- alignment/alignment-trimming - Pre-format trimming with column-mapping retention
- sequence-io/format-conversion - Convert sequence (non-alignment) formats

## References

- Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinf 29:2933-2935 (WUSS notation reference; see also the Infernal user guide).
- Larralde M et al. 2023. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinf 39:btad214.
- Cock PJA et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinf 25:1422-1423.
- Mistry J et al. 2021. Pfam: the protein families database in 2021. NAR 49:D412-D419.
- Steinegger M et al. 2019. HH-suite3 for fast remote homology detection and deep protein annotation. BMC Bioinf 20:473.
- Mirdita M et al. 2022. ColabFold: making protein folding accessible to all. Nat Methods 19:679-682.
- Blanchette M et al. 2004. Aligning multiple genomic sequences with the threaded blockset aligner. Genome Res 14:708-715.
