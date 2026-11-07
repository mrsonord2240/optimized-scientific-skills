## Subsample Reads (Deterministic, Pair-Consistent)

`samtools view -s SEED.FRAC` (same as `--subsample 0.FRAC --subsample-seed SEED`) -- integer is the hash seed; fractional is the keep fraction. The hash is on QNAME, so:
1. Mate consistency: read1 and read2 are kept or dropped together.
2. Reproducibility: same input file + same seed + same fraction returns the same reads. A bare `-s 0.1` is deterministic too (it means seed 0; two runs gave byte-identical output), but always write the seed explicitly: the seed used when none is given changed in samtools 1.24 (`--subsample 0.1` without `--subsample-seed` used seed 0 before 1.24; since 1.24 it derives the seed from a hash of the input header, `--subsample-seed auto`, and picks different reads: 1001 vs 1015 records on the test BAM).
3. **Sequential downsampling requires different seeds.** With the same seed the keep-sets are nested, so `-s 1.5` then `-s 1.25` keeps a nested 1/4 (25%) of the original, not 12.5%. Use different integer seeds for independent samples.

```bash
# 10% with seed 42 (always the same reads; pair-consistent)
samtools view -s 42.1 -b -o subset.bam input.bam

# Sequential cuts with INDEPENDENT seeds
samtools view -s 1.5 -b in.bam > half1.bam
samtools view -s 2.25 -b half1.bam > quarter.bam   # 12.5% of original
```

Coverage-matching to a target read count, or to another BAM's read count (tumor-normal: pull the tumor down to the normal). Hash-based, lands a few % off the target; a BAM already at or under the target is copied unchanged (`scripts/match_read_count.sh`):
```bash
bash scripts/match_read_count.sh input.bam matched.bam --target 10000000
bash scripts/match_read_count.sh tumor.bam tumor_matched.bam --like normal.bam
```

The script's guard matters: a fraction of 1 or more spliced into `-s` (`1.084419`) silently keeps only 8% of the reads.
