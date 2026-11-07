# DNA BAM, no UMIs: samtools markdup

```bash
ASSAY=wgs bash scripts/markdup.sh in.bam out.bam [threads]
```

`ASSAY` is required: `wgs wes capture somatic chipseq cutrun atac adna pacbio-amplicon`. The script refuses `rnaseq scrna umi amplicon longread 16s its`, and also refuses a BAM that looks spliced (STAR/HISAT2 in `@PG`, or N in CIGARs); `ALLOW_SPLICED=1` overrides that for a DNA BAM realigned by STAR. It collates, fixmates, sorts, marks, checks that no records were lost, indexes, and writes `out.markdup_stats.txt` beside the output.

Set the optical distance for the platform with `OPTICAL_DIST` (default 2500). Samtools' own default of 0 turns optical detection off.

| Platform | `OPTICAL_DIST` |
|----------|----------------|
| HiSeq 3000/4000/X, NovaSeq 6000, NovaSeq X, NextSeq 1000/2000 (patterned) | 2500 |
| HiSeq 2000/2500, MiSeq, NextSeq 500/550 | 100 |
| Element AVITI, MGI/DNBseq | different read-name format; supply a custom read-name regex |

- Read groups are respected (`--use-read-groups`); without it, independent libraries with the same coordinates are over-marked. For several lanes of one library, Picard keys on the LB tag instead (`routes/other-markers.md`).
- Re-marking an already marked BAM: add `-c` to the `samtools markdup` line in the script, or old flags survive.
- Optical (SQ) versus library (LB) duplicates: `samtools view -f 1024 out.bam | grep -o 'dt:Z:[A-Z][A-Z]' | sort | uniq -c`.
- Optical detection needs Illumina-style read names (tile and x:y in the name). SRA-renamed reads give a "cannot decipher read name" warning and no `dt:Z:SQ` tags, whatever `OPTICAL_DIST` is.
- Over 50% flagged prints a warning: recheck `ASSAY`.

Done when `out.bam` and `out.bam.bai` exist, the script printed the flagged percentage, and the count in equals the count out.
