# CRAM and related work

```bash
samtools view -T reference.fa input.cram
```

For offline CRAM use, populate a reference cache with `seq_cache_populate.pl -root CACHE_DIR reference.fa` and set `REF_PATH=CACHE_DIR/%2s/%2s/%s`; the cache directory must exist.

| Need | Skill |
|---|---|
| CRAM and SAM/BAM details | bio-sam-bam-basics |
| FASTA or alignment indexing | bio-alignment-indexing |
| BAM-reference validation | bio-alignment-validation |
| Pileup | bio-pileup-generation |
| VCF consensus | variant-calling/consensus-sequences (upstream only; not on this shelf) |

Done when CRAM decodes with its matching reference or the appropriate sibling task is selected.
