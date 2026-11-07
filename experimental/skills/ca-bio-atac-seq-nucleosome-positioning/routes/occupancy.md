# Occupancy and nucleosome calls

```bash
nucleoatac run --bed regions.bed --bam sample.bam --fasta genome.fa --out nucleo --cores 8
```

Use matched assembly and regions at least 500 bp; source guidance suggests at least 30M nuclear reads. NucleoATAC 0.3.4 is unmaintained and has a source-reported uninitialized accumulator that can produce empty nucleosome calls even after successful exit. Inspect `.nucpos.bed.gz`, not just exit status; without a validated fix report calls unavailable. For +1 assignment select the first call downstream of each strand-aware protein-coding TSS; check interval direction on both strands. See `references/calling.md` for limitations and the source-reported fix; scPrinter is an untested alternative in this skill.

Done when occupancy and valid position files are checked and region, depth and version are reported.
