# Remove primers

```bash
FWD=GTGYCAGCMGCCGCGGTAA REV=GGACTACNVGGGTWTCTAAT bash scripts/remove_primers.sh raw_reads trimmed
```

The example is 16S V4 515F/806R; supply the primers that made your library (degenerate IUPAC bases allowed, both required). The script runs `cutadapt -g FWD -G REV --discard-untrimmed` on every `*_R1_001.fastq.gz` and its `_R2_` mate, writing the same names into `trimmed/`. Other suffixes: `PATTERN_R1=` and `PATTERN_R2=`. Pairs without the primer are dropped, since a primerless read is suspect.

- QIIME2 equivalent: `qiime cutadapt trim-paired --p-front-f FWD --p-front-r REV --p-discard-untrimmed`.
- Primers are synthetic and often degenerate. Left on, they read as sequencing error, shift the truncLen frame and inflate chimeras (symptom: a large read fraction removed as chimeric).
- Install: `conda install -c bioconda cutadapt itsxpress`.

Done when `trimmed/` holds a trimmed pair per sample; then `routes/dada2-16s.md`.
