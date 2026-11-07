## What Flagstat Does Not Reveal

A 99% flagstat mapping rate does NOT mean the data is usable. Common false-positive scenarios:

1. **Adapter readthrough**: short fragments (insert < 2 * read_length) sequence into adapter; aligners soft-clip the adapter portion and flag the read as MAPPED. `samtools stats` has no soft-clip field on 1.24 (a `grep "bases soft-clipped"` silently returns nothing), so count soft-clipped bases from the CIGAR of primary mapped reads (>5% is a rule of thumb for adapter contamination; local aligners and split reads also soft-clip):
   ```bash
   samtools view -F 2308 input.bam | awk -F'\t' '
       {n++; c = $6
        while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
            len = substr(c, 1, RLENGTH-1) + 0; op = substr(c, RLENGTH, 1); c = substr(c, RLENGTH+1)
            if (op ~ /[MIS=X]/) q += len; if (op == "S") s += len}}
       END {if (!n || !q) {print "no primary mapped reads: nothing to compute" > "/dev/stderr"; exit 1}
            printf "soft-clipped bases: %d of %d (%.2f%%)\n", s, q, s/q*100}'
   ```
   Denominator: query bases (M, I, S, =, X) of primary mapped records. Exits 1 with a message on an empty or unmapped BAM.
2. **Off-target enrichment** (capture/WES): `picard CollectHsMetrics` PCT_OFF_BAIT or PCT_SELECTED_BASES. Interval lists need the reference dictionary header (`picard BedToIntervalList I=targets.bed O=targets.interval_list SD=ref.dict`):
   ```bash
   picard CollectHsMetrics I=input.bam O=hs_metrics.txt R=ref.fa \
       BAIT_INTERVALS=baits.interval_list TARGET_INTERVALS=targets.interval_list
   ```
3. **Low-complexity pile-up**: telomere/centromere reads mass at MAPQ-0; counted as mapped but useless. Detect via MAPQ distribution.
4. **Cross-sample contamination**: `verifybamid2` estimates FREEMIX (> 1% degrades somatic calling; > 5% breaks germline calling are commonly cited cut-offs); `somalier` checks sample identity/relatedness (`extract` then `relate`, and has a `contamination` subcommand). Both need a whole-genome or exome BAM: on a small slice they report "No reads found in any of the regions" / too few markers.
   ```bash
   # --SVDPrefix must include the .dat suffix of the panel files (e.g. 1000g.phase3.10k.b38.vcf.gz.dat)
   verifybamid2 --SVDPrefix panel.vcf.gz.dat --Reference ref.fa --BamFile input.bam --Output sample_vb
   # FREEMIX is in sample_vb.selfSM

   somalier extract -s sites.vcf.gz -f ref.fa -d extracted/ input.bam    # FASTA must contain the sites' contigs
   somalier relate extracted/*.somalier
   ```
   (Flags checked against `--help` of verifybamid2 2.0.3 and somalier 0.3.5; not run end to end here because no whole-genome BAM was available.)
5. **Wrong reference build**: a BAM aligned to GRCh37 viewed against GRCh38 looks fine to flagstat but produces nonsense pileups. Compare `@SQ M5:` from BAM header with `samtools dict ref.fa` (if the header carries no M5, compare `@SQ` names and lengths with `ref.fa.fai`) -- see alignment-validation.

## Insert Size Caveats

`samtools stats` reports the IS section for every pair with both mates mapped and splits the pairs into `inward oriented`, `outward oriented` and `other orientation` counts (checked on a synthetic mate-pair library with the proper-pair flag set and unset: `insert size average` 2000.0, 100 outward pairs both times). So:
- Mate-pair libraries (RF orientation): IS is reported, with outward-oriented counts dominating. The pysam snippets and `qc_report.py` only look at properly paired reads, so they report nothing when the aligner leaves the proper-pair flag unset
- `qc_report.py` drops templates >= `MAX_INSERT` = 8000 from its mean and median; `samtools stats -i 8000` (the default) instead counts them at 8000, so the two means differ on long-insert libraries (3555.5 vs 2286 on a test BAM with four templates of 8000 bp or longer)
- ATAC-seq: bimodal/multimodal expected (nucleosome ladder ~50/~180/~370 bp). Unimodal suggests poor transposition.
- RNA-seq: TLEN includes intron span -- mean meaningless
- Bisulfite (PBAT): orientation reversed; samtools may not flag proper pair
