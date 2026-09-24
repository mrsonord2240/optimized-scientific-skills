## Context-Specific Tags

Beyond the standard fields, downstream tools depend on optional tags whose presence depends on aligner and assay. Inspect with `samtools view input.bam | head -1 | tr '\t' '\n'` or pysam `read.get_tag('XX')`.

| Tag | Set by | Meaning | Required by |
|-----|--------|---------|-------------|
| NM:i | bwa, minimap2, samtools calmd | Edit distance to reference | Edit-distance filters, e.g. `samtools view -e '[NM]<=2'` |
| MD:Z | bwa, samtools calmd | Mismatch positions (text) | Not needed by `bcftools mpileup` or mapDamage (output identical with MD/NM stripped, checked); `samtools calmd` regenerates it |
| MC:Z | samtools fixmate -m (bwa mem also writes it) | Mate CIGAR | samtools markdup |
| ms:i | samtools fixmate -m | Mate score (lowercase per SAMtags); minimap2's own `ms:i` is an unrelated DP score | samtools markdup |
| RG:Z | aligner from -R | Read group ID | GATK BQSR, MarkDuplicates LB lookup |
| SA:Z | All split-read aligners | Other alignments of the read: `rname,pos,strand,CIGAR,mapQ,NM;` records (pos 1-based, each ends with `;`) | Sniffles, Manta, cuteSV, GRIDSS, Delly |
| NH:i | STAR, HISAT2 | Number of reported hits | featureCounts multimapper handling (checked: subread 2.0.6 -- default excludes all NH>1 alignments, `-M` includes them, `-M --fraction` weights each by 1/NH), Salmon (not verified here) |
| HI:i | STAR | Hit index among NH (1-based by default; `--outSAMattrIHstart 0` for 0-based) | Not RSEM: checked RSEM 1.2.28's `rsem-calculate-expression --help` and `convert-sam-for-rsem --help`, neither mentions HI; RSEM groups a multi-mapped read's alignments by consecutive same-QNAME lines, not by an HI tag |
| XS:A | STAR (`--outSAMstrandField intronMotif`), HISAT2 | Strand inferred from splice motif | StringTie, Cufflinks |
| ts:A | minimap2 `-ax splice` | Transcript strand from splice motif | StringTie |
| CB:Z | Cell Ranger, STARsolo | Corrected cell barcode; checked STARsolo 2.7.11b (`--soloType CB_UMI_Simple` + whitelist writes `CB:Z` matching the whitelist entry). Cell Ranger not run (registration-gated download) | scRNA quantification |
| UB:Z | Cell Ranger, STARsolo | Corrected UMI; checked STARsolo 2.7.11b (`UB:Z` present alongside `CB:Z` in the same run). Cell Ranger not run (registration-gated download) | UMI-aware dedup |
| RX:Z | fgbio AnnotateBamWithUmis | Raw UMI (bulk) | fgbio GroupReadsByUmi |
| MI:Z | fgbio GroupReadsByUmi | Molecular identifier (UMI group) | CallMolecularConsensusReads, duplex calling |
| cs:Z | minimap2 --cs | Compact CIGAR-with-bases | paftools, SV tools |

A missing tag can make a tool refuse or quietly do less. `samtools markdup` on a file without fixmate's tags refuses (`no ms score tag. Please run samtools fixmate on file first.`, exit 1); check each tool's tag requirements instead of assuming.

## Provenance: @PG Chain

The `@PG` lines record every tool that touched the BAM, linked through `PP` (previous program) tags. This is the audit trail.

```bash
samtools view -H input.bam | grep '^@PG'
```
`samtools view -H` appends its own `@PG` line to the output; add `--no-PG` to see the file's chain unchanged.

A clean germline pipeline:
```
@PG ID:bwa PN:bwa VN:0.7.17
@PG ID:samtools PN:samtools VN:1.20 PP:bwa CL:samtools sort
@PG ID:samtools.1 PN:samtools VN:1.20 PP:samtools CL:samtools fixmate
@PG ID:samtools.2 PN:samtools VN:1.20 PP:samtools.1 CL:samtools markdup
```

A broken/missing chain (no PP, unknown tools, gaps) means the BAM cannot be reliably reproduced.
