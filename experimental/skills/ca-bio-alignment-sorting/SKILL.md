---
name: ca-bio-alignment-sorting
category: Data Analysis
description: Use when sorting, merging or checking the record order of BAM or CRAM files before indexing, variant calling or paired-end steps.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

# Alignment Sorting

Pick the one row that matches what you have. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| You have | Read |
|----------|------|
| A BAM or CRAM to sort (coordinate, name, tag, template-coordinate) | `routes/sort.md` |
| FASTQ and a bwa index, want a sorted indexed BAM | `routes/align-sort.md` |
| Several BAMs to combine | `routes/merge.md` |
| A BAM that may be mis-sorted or mislabelled | `routes/check-order.md` |
| A downstream tool and no idea which sort it needs, or mates-adjacent FASTQ extraction | `routes/required-order.md` |
| Python only | `routes/python.md` |
| Duplicate marking, viewing or converting, indexing, filtering | `routes/other-skills.md` |
| An error, or a slow or large sort | `routes/errors.md` |

Rules on every route:

- Only coordinate-sorted files can be indexed. `-n`, `-N` and `-t` output fail `samtools index`.
- `samtools sort -n` is natural order (`read2` before `read10`); `-N` is ASCII. Anything read by Picard, GATK or htsjdk needs `-N`.
- The `@HD SO:` header is a hint. Verify the records (`routes/check-order.md`).
- Keep the original until `samtools quickcheck -v` passes on the output.
- Use `-u` between piped samtools commands; set `pipefail` on any pipe from an aligner.

Tested with samtools 1.24, pysam 0.24.1, Picard 3.5.0, bwa 0.7.19.
