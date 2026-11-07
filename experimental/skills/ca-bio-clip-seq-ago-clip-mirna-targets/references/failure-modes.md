# Failure modes and diagnostics

## Standard AGO CLIP cannot assign a miRNA directly

- **Trigger:** An AGO HITS-CLIP or eCLIP peak set is used to request a direct
  per-miRNA target list.
- **Mechanism:** The bound target fragment is recovered without the miRNA's
  identity.
- **Symptom:** Seed matching assigns many candidate miRNAs to each peak.
- **Action:** Use a chimeric library for direct pairing, or clearly label the
  assignments as computational and restrict them to expressed miRNAs.

## Chimera recovery is sparse

- **Trigger:** An unenriched chimeric library is expected to provide deep
  coverage for every miRNA.
- **Mechanism:** Chimeric molecules are a minority of total reads and are
  distributed across many miRNAs.
- **Symptom:** Rare miRNAs have very few recovered interactions.
- **Action:** Use a targeted miR-eCLIP enrichment strategy or increase depth
  only after comparing protocol-matched pilot recovery, controls, and the
  primary source. The inherited numeric depth and recovery examples lack
  claim-level provenance and are not operational thresholds.

## Short miRNA segments are missed

- **Trigger:** The selected alignment mode is insensitive to 21-23 nt miRNA
  segments.
- **Symptom:** Few chimeras are recovered and an advertised short-read mode
  recovers more candidates.
- **Action:** Verify and use the installed tool's supported short-read alignment
  route; do not assume the source's Hyb flag syntax is current.

## Computational predictions dominate

- **Trigger:** TargetScan or another prediction database is treated as ground
  truth without CLIP evidence.
- **Symptom:** Thousands of targets are reported for one miRNA.
- **Action:** Separate prediction-only candidates from sites that overlap an
  AGO peak and from directly observed chimeras.

## Non-canonical pairs are lost

- **Trigger:** Only canonical seed matches are retained.
- **Mechanism:** Some direct pairs rely on 3'-supplementary or other
  non-canonical duplex features.
- **Symptom:** Direct chimeras are discarded because TargetScan has no site.
- **Action:** Preserve the direct chimera and analyze the full duplex with a
  separately validated method.

## The expression filter is missing or mismatched

- **Trigger:** Database miRNAs are assigned without matched small-RNA-seq.
- **Symptom:** Results are dominated by miRNAs absent from the cell type or
  present only at very low expression.
- **Action:** Use matched organism, tissue or cell type, and condition; state
  the threshold and retain the underlying expression values.

## HEAP is transferred outside its model

- **Trigger:** A Halo-Ago2 transgenic-mouse workflow is proposed for human
  tissue without an equivalent engineered system.
- **Action:** Use an appropriate human AGO CLIP protocol and label HEAP evidence
  as mouse-model evidence.

## Weak 6mer sites dominate

- **Trigger:** All 6mer matches are included in the main target set.
- **Symptom:** The weak-site class overwhelms 7mer and 8mer evidence.
- **Action:** Report 6mer sites separately and keep stronger canonical classes
  in the primary computational overlap.

## Strand or coordinate context is lost

- **Trigger:** BED or tabular transformations drop strand, transcript, or
  assembly identity.
- **Symptom:** Target sites overlap the wrong transcript orientation or cannot
  be reconciled with a prediction set.
- **Action:** Preserve strand throughout, record assembly and transcript
  versions, and verify coordinate conversions before intersection. TargetScan
  8 UTR-relative sites must pass through a matching spliced transcript map;
  never label the source coordinates as genomic BED.

## Hyb assignments vary across clean reruns

- **Trigger:** The same read id receives different RNA partners in repeated
  executions, even with pinned inputs and versions.
- **Mechanism:** Multi-hit selection or tie ordering can change the selected
  pair. In the pinned source, Perl hash iteration feeds a count-only tie sort,
  so process-random key order can also reassign collapsed read ids; a single
  `*_hybrids_ua.hyb` file does not prove repeat stability.
- **Action:** Use the shipped wrapper's fixed hash, locale, timezone, and
  single-thread controls plus its all-replicate consensus parser. Inspect
  `support.tsv`; exclude `unstable_assignment` and `missing_from_replicate`
  rows and report their counts rather than selecting a convenient partner.

## Targeted miR-eCLIP UMI length is ambiguous

- **Trigger:** The pinned targeted CWL's 9-nt prose or the called script's
  10-nt default is treated as the library's UMI contract.
- **Mechanism:** Those pinned interfaces disagree, and the script's explicit
  length argument is not integer-converted.
- **Action:** Obtain the length from the versioned protocol or library record,
  then use `extract_targeted_umi.py` with explicit length, library id, and
  protocol source. Omission fails closed; retain the emitted manifest.

## Cross-context comparisons disagree

- **Trigger:** Mouse versus human, tissue versus cell line, or different miRNA
  expression contexts are compared as if interchangeable.
- **Action:** Stratify by biological context and matched expression rather than
  forcing one merged target list.

## Probe enrichment underperforms

- **Trigger:** A miR-eCLIP enrichment library underperforms its declared,
  protocol-matched control or prospectively defined recovery target.
- **Action:** Verify probe identity, specificity, affinity, and library design;
  test multiple probes and retain an unenriched or input comparison.
