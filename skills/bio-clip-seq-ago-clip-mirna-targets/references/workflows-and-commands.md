# Workflows and command patterns

The current executable route below was validated against public Hyb commit
`028ab6371ce793ca5e86f475fce1f2cc6ad3c677`. Protocol-specific routes remain
separate; confirm the exact library layout before preprocessing.

## eCLIP-style preprocessing

For the Yeo total chimeric-eCLIP route at public commit
`75fe74e90e6e4ca670a5af76836d80db09bdbcb1`, the CWL `extract_umi` surface
extracts a ten-nucleotide UMI from read 1. It does not state that the first ten
bases of read 2 are also a UMI:

```bash
umi_tools extract --bc-pattern=NNNNNNNNNN \
    --stdin=R1.fq.gz --read2-in=R2.fq.gz \
    --stdout=R1.umi.fq.gz --read2-out=R2.umi.fq.gz

cutadapt -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 6 -m 18 -o R1.trim.fq.gz -p R2.trim.fq.gz \
    R1.umi.fq.gz R2.umi.fq.gz
```

Use the adapter, barcode, read orientation, and quality rules from the actual
library protocol. Do not assume these example sequences or a ten-nucleotide UMI
apply to every chimeric library. In particular, targeted miR-eCLIP in the same
source uses a different read structure and takes its UMI from the read-2
prefix. The pinned CWL prose says 9 nt, while the called script defaults to 10
nt and cannot safely accept an explicit length because its argument is not
converted to an integer. Neither value is an operational default. Obtain the
UMI length from the actual library protocol and run the local fail-closed route:

```bash
python scripts/extract_targeted_umi.py \
  --read1 R1.fastq.gz --read2 R2.fastq.gz \
  --umi-length "$PROTOCOL_DECLARED_UMI_LENGTH" \
  --library-id "$LIBRARY_ID" --protocol-source "$VERSIONED_PROTOCOL_OR_RECORD" \
  --output-fastq targeted.umi.r1.fastq --manifest targeted.umi.manifest.json
```

The manifest binds the declared length and protocol source to the input and
output hashes. Omission of the length or provenance fails closed. The route
does not infer a length from either conflicting pinned-source statement. Do not
run the generic paired command above on the targeted layout. Record one of
these states before execution:

| Declared layout | UMI route | Required evidence |
|---|---|---|
| Yeo total chimeric-eCLIP | 10 nt at read-1 5' end via `umi_tools extract` | Pinned CWL/README commit and adapter set |
| Yeo targeted miR-eCLIP | Protocol-declared length from the read-2 prefix via `extract_targeted_umi.py`, plus separate primer handling | Library id, versioned protocol or library record, declared length, target miRNA, primer, source commit, and read structure |
| Other CLEAR-CLIP/CLASH/eCLIP | Protocol-specific; no default | Published/library protocol and an executed bounded read-layout check |

## Current Hyb contract

Current Hyb does not accept an arbitrary FASTA pathname as `db`. Build and
validate a named database below `HYB_HOME/data/db` (including a Bowtie2 index),
retain the reference FASTA hashes, and invoke an explicit goal and id. The
combined reference must use disjoint typed identifiers containing
`_microRNA_` for mature miRNAs and `_mRNA_` for target transcripts; the parser
fails closed when fields 4/10 do not contain exactly one of each type:

```bash
hyb \
    in=R1.trim.fq.gz \
    id=sample \
    db=mirna_transcripts_GRCh38_GENCODEv49 \
    qc=none \
    align=bowtie2 \
    type=mim
```

The expected current output is
`sample_comp_mirna_transcripts_GRCh38_GENCODEv49_hybrids_ua.hyb`. The shipped
wrapper fixes the process environment to C locale, UTC, zero Python and Perl
hash seeds, disabled Perl key perturbation, and one thread for Hyb and common
numeric runtimes. These controls stabilize the pinned source's otherwise
process-random Perl hash traversal before its count-only tie sort; they do not
choose among discordant biological assignments. The wrapper runs at least two
clean replicates, verifies the expected non-empty output from each, and then
applies the all-run consensus policy. It does not build into a shared
`HYB_HOME`; database preparation is an explicit, separately reviewed stage.

## Candidate-chimera inspection

For an existing AGO CLIP BAM, the source proposes counting soft-clipped
records as a diagnostic candidate pool:

```bash
samtools view -h dedup.bam | awk '$6 ~ /S/' | wc -l
```

This is not itself a chimera call. Supplementary alignments, clipping causes,
read layout, and pipeline-specific junction logic must be accounted for.

## Schema, expression, and deterministic aggregation

Current `.hyb` rows have 16 tab fields: read id, sequence, folding energy,
segment-1 id and five coordinate/score fields, segment-2 id and five
coordinate/score fields, and the optional final attribute. RNA ids are fields
4 and 10 and can occur in either orientation. The shipped
`consensus_hyb.py` validates that schema, requires one `_microRNA_` and one
`_mRNA_` typed id, preserves both segments' coordinates/scores and orientation,
and uses exact field equality for expression filtering.

Expression input must name `mirna_id`, `expression_value`, `expression_unit`,
and `expression_source`. A numeric threshold is required with the table. The
output `sites.tsv` preserves those values; `targets.tsv` reports stable read-id
support and unique target-coordinate sites. `support.tsv` reports every
normalized assignment observed for each read, its supporting run indices and
count, and whether it passed the all-run rule. These are recovery measurements,
not binding-affinity estimates.

## Selected-miRNA enrichment pattern

For an enrichment library targeting one miRNA, filter the structured field:

```bash
awk -F '\t' 'NR==1 || $4=="MIMAT0000076_MirBase_miR-21_microRNA"' sites.tsv \
  > mir21_sites.tsv
```

Use a field-aware filter to avoid substring collisions once the real schema is
known, and verify the expected enriched miRNA from the library metadata.

## TargetScan overlap pattern

Official TargetScan 8 site positions are relative to transcript 3' UTRs. They
are not genomic BED. First normalize the desired release to a TSV with
`transcript_id`, one-based inclusive `utr_start_1`/`utr_end_1`,
`mirna_family`, and `site_type`. Supply a version-matched spliced UTR map with
assembly, annotation release, TargetScan release, strand, and zero-based
half-open exon arrays. Convert to BED12:

```bash
python scripts/targetscan_sites_to_bed12.py \
  --sites targetscan8.normalized.tsv --utr-map utr-map.tsv \
  --assembly GRCh38 --annotation-release GENCODEv49 \
  --targetscan-release 8.0 --output targetscan8.GRCh38.bed12 \
  --manifest targetscan8.GRCh38.manifest.json

bedtools intersect -split -s -wa -wb \
  -a chimera_sites.GRCh38.bed -b targetscan8.GRCh38.bed12 \
  > chimera_targetscan_supported.tsv
```

The converter fails on missing transcripts, release mismatches, out-of-range
sites, invalid exon arrays, or empty output. Verify chromosome naming before
intersection. Call a valid overlap
prediction-supported direct evidence rather than using the overlap alone as
functional validation.

## Expected inputs and outputs

The wrapper requires a non-empty trimmed FASTQ, a validated named Hyb database,
a safe run id, and a new output directory. Optional expression filtering uses
the exact schema above. It refuses existing output unless `--replace` is
explicit. Work occurs in an owned sibling staging directory; a failed Hyb or
parser stage returns nonzero and publishes nothing. A successful run atomically
publishes raw per-replicate Hyb outputs/logs, structured sites and targets,
per-read assignment support, excluded reads with reason codes, and a JSON
manifest with hashes, counts, and deterministic process controls.
