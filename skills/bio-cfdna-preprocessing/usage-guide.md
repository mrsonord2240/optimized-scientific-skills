# cfDNA Preprocessing - Usage Guide

Use this skill to choose and execute plasma cfDNA preprocessing without
destroying fragment-length signal or treating nucleosome-coincident molecules
as ordinary coordinate duplicates.

## Prerequisites

```bash
conda install -c bioconda fgbio bwa samtools
pip install pysam numpy
```

The bundled pipeline also needs a reference FASTA and an unmapped, UMI-bearing
BAM with a known read structure. Confirm versions and live CLI help before
running a real sample.

## Example requests

- "Build UMI consensus reads from my targeted cfDNA panel BAM."
- "Assess whether duplex consensus improves my validated low-VAF assay."
- "Explain the correct fgbio read structure for a UMI followed by a stem."
- "Do minimal preprocessing for sWGS tumor-fraction estimation."
- "Check whether my insert-size distribution suggests gDNA contamination."
- "My library has no UMIs; how should I handle duplicates quantitatively?"
- "Why is my adaptase library mode about 10 bp shorter than 167 bp?"

## Start here

1. Read the core workflow in [SKILL.md](SKILL.md).
2. Choose chemistry and consensus depth with [method selection and fragment
   QC](references/method-selection-and-fragment-qc.md).
3. Use the exact [fgbio consensus workflow](references/fgbio-consensus-workflow.md)
   when assembling CLI commands.
4. Import [`scripts/preprocess_cfdna.py`](scripts/preprocess_cfdna.py) when a
   reusable Python wrapper or insert-size summary is appropriate.
5. Cite the claims using the [scientific references](references/scientific-references.md).

The essential gates are: map every molecular UMI segment to a tag, use paired
grouping for duplex data, call consensus permissively and filter afterward,
query-group before template-aware filtering, coordinate-sort/index only after
filtering, never use naive coordinate deduplication as a quantitative no-UMI
cfDNA strategy, and do not size-select before reporting fragmentomics.

The Python wrapper accepts paths containing spaces and shell metacharacters as
literal filenames. It rejects missing inputs, unsafe input/output aliasing,
invalid UMI tag layouts, and thread counts outside 1-256 before launching tools.

## Focused regression suite

Point `CFDNA_FIXTURE_ROOT` at the prepared bounded fixture directory, then run:

```bash
python -m unittest discover -s tests -v
```

The live suite covers simplex and reciprocal-duplex extraction, ordinary and
space/metacharacter paths, final coordinate order/indexing, QC flag/boundary
behavior, empty input, and ten consecutive complete wrapper calls.
