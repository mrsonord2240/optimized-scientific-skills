# UMI consensus

```bash
python scripts/cfdna_cli.py consensus raw.unmapped.bam final.bam reference.fa
```

Run from the skill directory with paired UMI-bearing unmapped BAM, indexed reference and writable intermediate space. Defaults are simplex, eight threads, two `6M11S+T` read structures and tags ZA/ZB; change these for the actual library using `--read-structures`, `--tags`, and `--threads`. Add `--duplex` only for a suitable duplex design. The wrapper calls fgbio, bwa and samtools and requires numpy and pysam. It calls consensus permissively, re-aligns, then filters; simplex final minimum is 2 reads, duplex is 2 total with 1 on each strand. Do not assume that this fixed filter retains singletons. Check installed tool help and verify assay-specific choices before use; intermediate files are written beside the output and existing names may be overwritten. The detailed [manual workflow](../references/fgbio-consensus-workflow.md) covers alternative policies and flags.

The CLI delegates to the retained source implementation `scripts/preprocess_cfdna.py`.

**Done when:** coordinate-sorted indexed output and the chemistry, read structure, grouping, thresholds, family counts, tool versions and reference are recorded.
