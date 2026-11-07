# WASP mapping-bias filtering

```bash
bowtie2 -x hg38 -1 R1.fq -2 R2.fq -S step1.sam
```

This is the first step of the complete multi-command [WASP recipe](../references/methods-and-commands.md), not a finished ASB test. Adapt the index and input paths; generate sample-specific heterozygous SNP tables, align and sort, swap and remap paired reads using the same alignment settings, filter inconsistent remappings, then merge, sort and index. Verify filenames and flags against the installed WASP release. Record input, direct-keep and remap-keep counts and verify paired output; do not assume a fixed loss fraction. Follow with [BaalChIP](baalchip.md) or another suitable model.

**Done when:** a sorted indexed bias-filtered BAM and stage counts are retained, with reference, sample variants and tool versions documented.
