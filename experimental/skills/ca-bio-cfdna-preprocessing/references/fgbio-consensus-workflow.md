# fgbio consensus workflow

Use this recipe for targeted cfDNA libraries carrying inline single-strand or
duplex UMIs. Verify every flag against the installed fgbio version before
adapting it to a production pipeline.

## Required order

Consensus needs alignment coordinates for grouping, but the called consensus
sequence is emitted unmapped. The valid order is therefore:

`extract UMI -> query-group -> FASTQ -> align/zip tags -> template-sort -> group -> call consensus -> re-align/zip tags -> query-sort -> filter -> coordinate-sort/index`

```bash
# 1. Extract each inline UMI segment into a molecular-index tag and also
#    combine them in RX. M=UMI, S=skip/stem, T=template, +=all remaining.
fgbio ExtractUmisFromBam --input raw.unmapped.bam --output with_umis.bam \
    --read-structure 6M11S+T 6M11S+T \
    --molecular-index-tags ZA ZB --single-tag RX

# 2. Query-group the uBAM, emit interleaved paired FASTQ, align, then zipper
#    the uBAM metadata/UMI tags back onto mapped records. Paths are quoted.
samtools sort -n --threads 8 -o with_umis.queryname.bam with_umis.bam
samtools fastq -T RX,ZA,ZB with_umis.queryname.bam \
  | bwa mem -C -t 8 -K 150000000 -Y -p "reference.fa" - \
  | fgbio ZipperBams --unmapped "with_umis.queryname.bam" \
      --ref "reference.fa" --output aligned.unsorted.bam
samtools sort --template-coordinate --threads 8 \
    -o aligned.template-coordinate.bam aligned.unsorted.bam

# 3. Group. Use paired for duplex; use adjacency for single-strand UMIs.
fgbio GroupReadsByUmi --input aligned.template-coordinate.bam --output grouped.bam \
    --strategy paired --edits 1

# 4a. SIMPLEX: --min-reads takes one value.
fgbio CallMolecularConsensusReads --input grouped.bam \
    --output consensus.unmapped.bam --min-reads 1

# 4b. DUPLEX: keep the caller pre-filter permissive and filter later.
fgbio CallDuplexConsensusReads --input grouped.bam \
    --output consensus.unmapped.bam --min-reads 1

# 5. Re-align, transfer consensus/UMI tags, and queryname-sort for the
#    template-aware consensus filter.
samtools fastq consensus.unmapped.bam | bwa mem -t 8 -Y -p reference.fa - \
  | fgbio ZipperBams --unmapped consensus.unmapped.bam --ref reference.fa \
      --tags-to-reverse Consensus --tags-to-revcomp Consensus \
      --output consensus.mapped.unsorted.bam
samtools sort -n --threads 8 -o consensus.mapped.queryname.bam \
    consensus.mapped.unsorted.bam

# 6. Filter query-grouped records, then coordinate-sort and index. The
#    min-reads values are total, strand1, strand2.
fgbio FilterConsensusReads --input consensus.mapped.queryname.bam \
    --output filtered.queryname.bam \
    --ref reference.fa --min-reads 2 1 1 \
    --max-read-error-rate 0.025 --max-base-error-rate 0.1 \
    --min-base-quality 40 --reverse-per-base-tags
samtools sort --threads 8 -o filtered.bam filtered.queryname.bam
samtools index --threads 8 filtered.bam
```

The simplex and duplex caller commands are alternatives; do not run both on
the same grouped BAM unless intentionally producing separate outputs.

## Flag semantics to preserve

- `--strategy paired` is mandatory for duplex because `adjacency` cannot
  reconstruct A/B strand pairing when the UMI appears in opposite order.
- `CallDuplexConsensusReads --min-reads` is treated here as a permissive
  pre-filter. Use `1` and apply the specificity threshold later.
- `FilterConsensusReads --min-reads 2 1 1` requires two reads total with both
  strands represented. A policy such as `1 1 0` permits a missing strand.
- If strand thresholds differ, put the more stringent value first; for example
  use `6 3 0`, not `6 0 3`.
- `--reverse-per-base-tags` orients per-base depth/error tags consistently with
  genomic orientation after alignment.
- Preserve the `S` token in read structures with a non-template stem. Omitting
  it allows stem sequence to bleed into the template.
- Supply exactly one unique two-character SAM tag for every `M` segment;
  `--single-tag RX` combines UMIs but does not replace the per-segment mapping.
- `ZipperBams` requires the mapped and unmapped inputs to have identical
  queryname grouping and order. Query-group before FASTQ conversion.
- `FilterConsensusReads` likewise needs query-grouped templates. Coordinate
  sort and index only after filtering.

## Execution-specific failure modes

| Symptom | Likely cause | Required response |
|---|---|---|
| Consensus BAM has invalid or missing coordinates | Consensus output was not re-aligned | Re-align, transfer tags with `ZipperBams`, then filter |
| UMI/stem bases corrupt alignment | Read structure omitted `S` | Encode the actual stem, for example `6M11S+T` |
| Duplex run yields no two-strand families | Data grouped with `adjacency` | Regroup with `--strategy paired` |
| Molecular recovery collapses before filtering | Duplex caller threshold set too high | Call with `--min-reads 1`; filter afterward |
| Strand specification is rejected or misapplied | Strand values are in the wrong order | Confirm installed help and place the stricter strand first |
| Consensus tags disappear after mapping | New alignment was not zipped to unmapped consensus | Include `ZipperBams` before sorting the mapped consensus |
| Extraction reports molecular segments but zero tags | `--single-tag RX` was treated as a segment mapping | Add one `--molecular-index-tags` value per `M` segment and keep `RX` only as the combined tag |
| BWA reads a few long garbage sequences or emits invalid SAM | A BAM pathname was supplied where FASTQ was required | Query-group the uBAM, stream `samtools fastq` to `bwa mem -p`, and zipper metadata back |
| Consensus filtering rejects sort order | Coordinate-sorted input was supplied | Queryname-sort before filtering; coordinate-sort/index the filtered output |

The bundled [`../scripts/preprocess_cfdna.py`](../scripts/preprocess_cfdna.py)
encodes the same sequence as importable Python. Its external command surfaces
still require live-tool verification and representative sample testing.
