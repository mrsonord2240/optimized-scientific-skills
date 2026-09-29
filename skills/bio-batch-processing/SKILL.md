---
name: bio-batch-processing
description: Process many sequence files in batch (count, merge, split, convert, summarize) with memory-safe streaming and on-disk indexing using Biopython, pysam, or pyfastx. Use when iterating over a directory of FASTA/FASTQ files, merging or splitting datasets, building random access across many or huge files, or automating per-file operations without exhausting RAM.
tool_type: python
primary_tool: Bio.SeqIO
license: MIT
author: GPTomics
---

# Batch Processing

Process directories of sequence files without materializing the complete dataset in memory.

## Version compatibility

The source examples target Biopython 1.83+, with pysam 0.22+ and pyfastx 2.0+
as optional readers. Before using a pattern, check the installed package:

```bash
pip show biopython
pip show pysam pyfastx
```

If an import, attribute, or signature differs, inspect the installed API with
`help(module.function)` and adapt the call. Do not retry an incompatible example
unchanged.

## Core rules

- Stream with `SeqIO.parse()`; do not wrap a large parse in `list()`.
- Count with `sum(1 for _ in records)` rather than `len(list(records))`.
- Re-create a parse generator for each pass because it is one-shot.
- Sort discovered paths with the shared stable relative-path key before every
  multi-file operation; filesystem traversal order is not reproducible.
- Use `SeqIO.index_db()` instead of `to_dict()` for persistent random access
  across many or large files.
- Use `pysam.FastxFile` for thin linear FASTQ iteration at very large scale;
  use pyfastx when indexed reuse of plain or gzipped FASTA/FASTQ is required.
- Decide how to handle duplicate record ids before merging. A streamed FASTA may
  contain repeated ids, but `index_db()` and `to_dict()` reject duplicate keys.
- Treat format conversion as potentially lossy. GenBank-to-FASTA drops features,
  annotations, and qualifiers.

## Workflow

1. Identify input formats, compression, recursive-search needs, expected record
   count, and whether ids are unique across files.
2. Decide between one-pass streaming and persistent random access. Read
   [Reader selection and indexing](references/readers-and-indexing.md) for the
   reader matrix, compression constraints, and quality-encoding caveat.
3. Select the operation and use the routed recipe in
   [Batch operation recipes](references/operations.md).
4. Keep transforms generator-based. A bounded split batch may hold only the
   configured chunk, never the full input.
5. Close indexes and output handles, and report output paths, counts, and any
   rejected duplicate ids or lossy conversions.
6. For a prefix split, use a new explicit output directory and retain its
   completion manifest. Treat sequence ids as untrusted path input.

## Routing

| Need | Read or run |
|---|---|
| Choose SeqIO, `index_db`, pysam, or pyfastx | [Reader selection and indexing](references/readers-and-indexing.md) |
| Build/reuse a pyfastx gzip index and check random access | `python scripts/pyfastx_index.py INPUT.fa.gz --record-id ID` |
| Count files, merge, or tag source filenames | [Count and merge](references/operations.md#count-and-merge) |
| Split by record count or id prefix | [Split files](references/operations.md#split-files) |
| Convert formats | [Batch conversion](references/operations.md#batch-conversion) |
| Parallelize per-file work | [Parallel processing](references/operations.md#parallel-processing) |
| Produce per-file CSV summaries | [Summary statistics](references/operations.md#summary-statistics) |
| Run the standalone count/split/index demonstration | `python scripts/batch_process.py` |
| Run focused safety and portability regressions | `python -B tests/test_batch_process.py` |

## Minimal streaming pattern

```python
from pathlib import Path
from Bio import SeqIO

def stable_key(path):
    relative = path.relative_to("data").as_posix()
    return relative.casefold(), relative

for fasta_file in sorted(Path("data").glob("*.fasta"), key=stable_key):
    with fasta_file.open("r", encoding="utf-8") as handle:
        count = sum(1 for _ in SeqIO.parse(handle, "fasta"))
    print(f"{fasta_file.name}: {count} sequences")
```

Use `Path.rglob()` instead of `glob()` only when recursive discovery is intended.

## Failure checks

- A killed process or `MemoryError` usually means a parse was materialized.
- Slow plain iteration over tens of millions of reads may justify a thinner
  pysam reader.
- `ValueError: Duplicate key` means ids collide across the indexed inputs.
- An empty second loop usually means a parse generator was already exhausted.
- A rejected prefix split means an id is not a portable filename, two prefixes
  collide case-insensitively, or an output already exists. Inspect the partial
  manifest and start a reviewed retry in a new directory; do not delete or
  overwrite outputs automatically.
- `records_per_file must be a positive integer` means no records were consumed
  and no chunk was created; correct the requested chunk size and retry.
- Plain gzip is not seekable by `SeqIO.index_db()`; use BGZF for indexed
  Biopython access, or choose a reader that explicitly indexes gzip.
- Missing annotations after conversion are data loss, not a parser failure.

## Related skills

- `read-sequences` for per-file parse and index semantics
- `filter-sequences` for streamed per-record filtering
- `sequence-statistics` for N50 and length distributions
- `format-conversion` for conversion-specific data-loss traps
- `compressed-files` for BGZF versus plain gzip
- `paired-end-fastq` for synchronized mate processing
- `database-access/entrez-fetch` for NCBI batch downloads

## Provenance

Normalized from `GPTomics/bioSkills` at commit
`d91ed3d563019e649dc854c56ccd62551359488a`, path
`sequence-io/batch-processing`. Original author: GPTomics. License: MIT.
