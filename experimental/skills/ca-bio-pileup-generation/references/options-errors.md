## Pileup Options and Defaults

| Option | Description | Default (samtools / bcftools mpileup) | Common pitfall |
|--------|-------------|---------------------------------------|----------------|
| `-f FILE` | Reference FASTA | required by bcftools (`--no-reference` to skip) | Triggers BAQ ON by default; a missing file or contig gives reference base `N` |
| `-r REGION` | Restrict to region | | |
| `-l FILE` | BED file of regions | | |
| `-q INT` | Min mapping quality | 0 / 0 | Aligner-dependent semantics |
| `-Q INT` | Min base quality (after BAQ) | **13 / 1** | Silently drops low-quality bases, and low-quality `*` deletion slots; `-Q 0` with default overlap detection has subtle behavior |
| `--ff FLAGS` (bcftools `--ns`) | Skip reads with any of these flags | **UNMAP,SECONDARY,QCFAIL,DUP** | Pre-flagged duplicates vanish from depth; `--ff 0` keeps all mapped reads |
| `--rf FLAGS` (bcftools `--lu`) | Keep only reads with any of these flags set | none | bcftools `--nu` is not this: it skips a read missing any listed bit, i.e. keeps only reads with **all** bits |
| `-d INT` | Max depth per file | **8000 / 250** | Silently truncates; `-d 0` = no cap |
| `-B` | Disable BAQ | BAQ on with `-f` | Often correct for long reads, SV, viral, consensus; cannot combine with `-E` |
| `-A` | Count anomalous pairs (paired, not proper) | dropped | Required for amplicon (reads are by design not properly paired) |
| `-a` / `-aa` | `-a`: every position of each contig that has reads, incl. zero coverage; `-aa`: also contigs with no reads | off | Required for ARTIC SARS-CoV-2 consensus generation; `-a` = `-aa` for a single-contig reference |
| `-x` (`--disable-overlap-removal`; bcftools: `--ignore-overlaps`) | Disable mate-overlap correction | overlaps counted once | Rarely correct |
| `--max-BQ INT` (`bcftools mpileup` only) | Cap baseQ/BAQ | 60 | Not a `samtools mpileup` option; useful for ONT/HiFi (Q values inflated) |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| `[E::faidx_adjust_position] The sequence "X" was not found`, **exit status 0**, output rows with reference `N` | `-f` FASTA lacks the BAM's contig (`MN908947.3` vs `MT192765.1`, `chr1` vs `1`) | Compare `samtools view -H in.bam \| grep '^@SQ'` with `cut -f1 ref.fa.fai` before running; do not trust the exit status |
| Reference column all `N`, no error | `samtools mpileup` run without `-f` | Add `-f reference.fa` (`bcftools mpileup` refuses: "requires the --fasta-ref option") |
| `[E::mpileup] fail to parse region '22:3000-3005'` | Contig name mismatch (`chr22` vs `22`) | Use the name from the `@SQ` header |
| `[E::idx_find_and_load] Could not retrieve index file for 'in.bam'` | `-r` (or `bam.pileup()`) on an unindexed BAM | `samtools index in.bam` |
| Empty output | Region has no reads | `samtools view in.bam chr1:1000-2000 \| head`; check `ls reference.fa.fai` |
| Out of memory | High coverage region | Cap with `-d`, set well above the expected coverage (reads above the cap are dropped); or restrict with `-l targets.bed`; process contigs in parallel |
| `Failed to read from standard input: unknown file type` | Text `samtools mpileup` piped into `bcftools call` | Use `bcftools mpileup` |

Text pileup is for debugging and small regions; BCF is far cheaper to process at scale.
