### File Integrity
```bash
# Fast: header + EOF block check (misses mid-file truncation, invalid CIGAR)
samtools quickcheck -v in.bam || echo "QUICKCHECK FAILED"
samtools quickcheck -v *.bam > bad_bams.fofn   # one fail-line per bad file

# Slow but thorough: structural validation (R= enables the NM check; use samtools calmd to check MD)
picard ValidateSamFile I=in.bam MODE=SUMMARY R=ref.fa
```

What each check caught in 21 planted defects (synthetic BAMs derived from a real human PE BAM; samtools 1.24, Picard 3.5.0), and no flags on the 2 valid synthetic controls:

| Check | Caught | Misses |
|-------|--------|--------|
| `samtools quickcheck` | 3/21: missing EOF block, truncated file, header without `@SQ` | everything inside the file: bit flips, dropped blocks, bad records |
| `samtools view -c` (decodes every record) | 3/21: truncated file, bit-flipped block, CIGAR/SEQ length mismatch | dropped blocks and every semantic defect |
| `picard ValidateSamFile` | 18/21: mate fields (`MATE_NOT_FOUND`, `MISMATCH_MATE_ALIGNMENT_START`, `MISMATCH_FLAG_MATE_*`), unsorted (`RECORD_OUT_OF_ORDER`), read group missing from header, MAPQ on unmapped reads, CIGAR off the contig end, CIGAR/SEQ length, bit flips, a dropped mid-file block (as orphaned mates) | 0x40 and 0x80 both set, TLEN inconsistent with mate positions, header-only file with 0 records |
| CI one-liner below | 6/21: the quickcheck and decode catches plus 0-record files | all record-level defects |

No tool here sees a whole block removed from the middle unless mates are orphaned; for transfers compare `md5sum` or the read count with the source.

Picard noise on valid files (observed): `MATE_NOT_FOUND` on region-extracted BAMs (mates fall outside the slice), `MISSING_TAG_NM` warnings on BAMs without NM tags (STAR output), `RECORD_OUT_OF_ORDER` on a samtools name-sorted BAM (Picard likely orders query names differently from samtools). Long-read BAMs can also report missing `@RG PL` / other required header tags and NM-convention errors. Read the SUMMARY before trusting a non-zero count.

`IGNORE=<TYPE>` removes a check. On the planted files `IGNORE=INVALID_MAPPING_QUALITY IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND` turned 40, 2 and 2820 genuine errors into "No errors found". Ignore a type only after tracing it to a known harmless source.

CI-safe one-liner (`MIN_READS` = smallest read count you accept; default 1):
```bash
MIN_READS=${MIN_READS:-1}
test -s in.bam \
  && samtools quickcheck -v in.bam \
  && [ "$(samtools view -c -F 2304 in.bam)" -ge "$MIN_READS" ] \
  || { echo "BAM failed integrity"; exit 1; }
```

### Sequence Dictionary Cross-Validation (SN / LN / M5)
```bash
# Compare contig names and lengths, and M5 wherever the BAM header carries it. Exit 1 on any difference.
awk -F'\t' '
  function parse(   i,k,v) { sn=ln=m5=""; for (i=2;i<=NF;i++) { k=substr($i,1,2); v=substr($i,4); if (k=="SN") sn=v; else if (k=="LN") ln=v; else if (k=="M5") m5=v } }
  FNR==NR { if ($1=="@SQ") { parse(); rlen[sn]=ln; rm5[sn]=m5 } ; next }
  $1=="@SQ" { parse(); if (m5!="") nm5++
      if (!(sn in rlen))      { print "NOT IN REFERENCE:", sn; bad=1 }
      else if (rlen[sn]!=ln)  { print "LENGTH DIFFERS:", sn, ln, "vs", rlen[sn]; bad=1 }
      else if (m5!="" && m5!=rm5[sn]) { print "M5 DIFFERS:", sn; bad=1 } }
  END { if (!nm5) print "no M5 in BAM header: only names and lengths were compared"; exit bad }
' <(samtools dict ref.fa) <(samtools view -H in.bam)
```

Names and lengths (SN, LN) catch renaming (UCSC `chr1` vs Ensembl `1`), a different contig set and different lengths. M5 catches a different sequence under the same name and length, and a contig-to-sequence mix-up. Most aligners (bwa, minimap2, STAR) write no M5, so a BAM without it can only be checked on SN/LN; the snippet says so instead of failing. Concrete failure modes: GRCh38 vs GRCh38.p13 vs GRCh38_no_alt (alt contigs differ); soft-masked vs unmasked (M5 matches -- case is normalized to uppercase before hashing, so lowercase soft-masking is invisible to M5). Hard-masking is different: it replaces bases with `N`, changing the sequence, so its M5 does NOT match the unmasked reference. Where M5 is present it is the only definitive identity check.
