## HGVS Splicing Nomenclature

Following den Dunnen 2016 *Hum Mutat*:

| Notation | Meaning |
|----------|---------|
| `c.123+1G>A` | +1 of intron downstream of exon ending at cDNA position 123 (canonical 5'ss G) |
| `c.123+5G>A` | +5 position of donor (consensus region) |
| `c.124-1G>A` | -1 of acceptor (canonical AG) |
| `c.124-3T>G` | -3 of acceptor (Py-tract / BPS region) |
| `c.124-50A>G` | Deep-intronic; may activate cryptic site |
| `r.123_456del` | RNA-level deletion (predicted exon skipping) |
| `r.spl?` | Unknown splice consequence |
| `r.0?` | No detectable RNA |
| `p.0?` | Unknown protein consequence |
| `p.(=)` | No predicted protein change (silent) |

Validation tools: VariantValidator (Freeman 2018 *Hum Mutat*), Mutalyzer 2 (Lefter et al 2021 *Bioinformatics* 37:2811-2817). Check REF against the reference FASTA before scoring; a wrong REF is skipped without failing the run.

## Extended-Window Scoring for Deep-Intronic Variants

`-D` is the distance from the variant to the gained/lost site, not the distance to the nearest canonical site. A pseudoexon-creating variant usually scores **at the variant itself** (it creates a donor or acceptor there); a wider `-D` adds the *other* pseudoexon boundary. Measured on GLA c.639+919G>A (919 nt from exon 4): `-D 50` DS_DG 0.30 at the variant; `-D 500` and `-D 2000` add DS_AG 0.22 at +53 (the far end), delta_max unchanged at 0.30. Pangolin `-d 500 -s 0.2 -m False` lists both ends (+0.23 at -3 and +0.23 at +53).

| Window | Tradeoff |
|--------|----------|
| -D 50 (default) | Fast; captures canonical-site disruption and gains created at the variant |
| -D 500 | Adds far-end boundaries and distant sites (e.g. DMD c.31+1G>A donor gain 0.08 at -10 with -D 50 becomes 0.40 at +69) |
| -D 2000 | Maximum sensitivity; more distant, weaker sites |

For unsolved cases run `-D 500` on every intronic candidate that scores below 0.20 (including 0.00), and report DS **and DP of both ends** of any candidate pseudoexon:

```bash
spliceai -I candidates.vcf -O output_D500.vcf -R genome.fa -A grch38 -D 500 -M 0
```

Pseudoexon creation in deep introns explains a substantial fraction of unsolved Mendelian disease alleles in current cohorts (estimates 5-15% across studies; specific quantitative range will vary by cohort and panel — verify against current literature). Disease examples: CFTR 3849+10kbC>T, USH2A c.7595-2144A>G, CEP290 c.2991+1655A>G (LCA10), GLA c.639+919G>A (Fabry). Extended context is not a rescue: in the checked GRCh38 panel CFTR 3849+10kbC>T remained SpliceAI 0.16 at D50/D500/D2000, while Pangolin was +0.33 and MMSplice had no row. Treat this 0.10-0.20 pattern as discordant/inconclusive and seek RNA evidence. CI-SpliceAI (see above) is a second opinion, not a rescue.

## Concordance Across Predictors

`build_concordance()` in `scripts/splice_parsers.py` parses the three tools' real outputs, joins them on a normalised `chrom:pos:ref>alt` key (`chr` prefix dropped), keeps every input variant (outer join: a variant a tool skipped shows NaN and lowers `n_scored`) and never drops a tool's missing output silently:

```python
from scripts.splice_parsers import build_concordance
T = build_concordance('input.vcf', spliceai_vcf='output.vcf', pangolin_vcf='pangolin_output.vcf',
                      mmsplice_csv='mmsplice_predictions.csv')
print(T[['spliceai_delta', 'pangolin_score', 'delta_logit_psi', 'n_scored', 'n_above', 'concordance']])
```

Thresholds per tool: SpliceAI delta >= 0.2, Pangolin |score| >= 0.2, MMSplice |delta_logit_psi| >= 1.0 (Skill conventions for Pangolin/MMSplice; only the SpliceAI one is ClinGen). Multi-allelic records must be split first (`bcftools norm -m-`): Pangolin scores ALT[0] only.

| Concordance label | Meaning | Action |
|-------------------|---------|--------|
| `all_predict_disruption` | every tool that scored is above threshold | Report `n_above/n_scored` and score margins; candidate supporting evidence only after transcript, population and RNA review |
| `majority_predict_disruption` | more than half above (e.g. 2/3) | Report `n_above/n_scored`; do not upgrade evidence strength from agreement alone |
| `discordant` | some, but not most, above | Report inconclusive; flag for RNA validation |
| `none_predict_disruption` | all below | BP4 (supporting) only if SpliceAI <= 0.10 |
| `insufficient_tools` | fewer than 2 tools scored | Find out why (skipped records) before interpreting |

Discordance is the most informative pattern — variants where one model sees impact and others don't are high priority for RNA validation. Always show `n_scored`, `n_above`, and the underlying scores: on the test panel GLA c.639+919G>A is `all_predict_disruption` from only two tools because MMSplice returned nothing. Concordant false positives occur, so never infer pathogenicity or RNA outcome from this label.

## Branchpoint Variant Detection

All current tools are **weak at branchpoint variants** because the BPS motif (yUnAy) has low information content; SpliceAI captures only some. Published branchpoint-specific methods include BPP, LaBranchoR, SVM-BPfinder and BPHunter. **None was run here.** The BPHunter standalone reference URL redirected to a GitHub 404 when checked on 2026-09-20; do not present it as an available local workflow. For a borderline SpliceAI delta (0.1-0.3) in the BPS region (-18 to -40 from 3'ss), treat predictors as uninformative and require RNA validation; an available web service may be used only as a supplementary research check.

## Splice-Switching ASO Design

**Goal:** Plan antisense oligonucleotides to modulate splicing therapeutically (e.g. SMA ISS-N1, DMD exon skipping). This is a **checklist for discussion with a design platform, not an executable pipeline**; no design tool is provided or was run here, and SpliceAI cannot model oligonucleotide occlusion (on N-masked sequence it returns `.` for every score).

1. Identify the target ESE/ESS/ISE/ISS region (splice-site strength with MaxEntScan; motif databases)
2. Design candidate 18-22 nt ASOs spanning the regulatory element
3. Filter for RNA accessibility (avoid stable hairpins) with RNAfold
4. Transcriptome-wide sequence search for off-target binding (<=16/20 nt match to any non-target transcript)
5. Avoid TLR9 immunostimulatory CpG motifs
6. Score the *variant/sequence change* you aim to correct or mimic with SpliceAI, never the oligo

Chemistry choices: 2'-MOE-PS (nusinersen-like, CNS, intrathecal); PMO (DMD ASOs, systemic IV); GalNAc-conjugated (hepatic targeting). Approved precedents: **nusinersen** (SMA ISS-N1 occlusion, exon 7 inclusion); **risdiplam** (small-molecule SMN2 splicing modulator); **eteplirsen/golodirsen/casimersen/viltolarsen** (DMD exon skipping). Design references: Hua 2008 *AJHG*; Roberts et al 2023 *Nat Rev Drug Discov* 22:917 (DMD therapeutic approaches).

## Per-Tool Failure Modes

### SpliceAI: Effect Site Farther Than -D

**Trigger:** The gained/lost splice site lies more than `-D` nt from the variant (deep-intronic pseudoexon whose creating change is far from its boundaries, or a distant cryptic site).

**Mechanism:** `-D` limits the variant-to-site distance; nothing outside is scored.

**Symptom:** Known pathogenic deep-intronic variant scores low (<0.2) or 0.00; no pseudoexon detected.

**Fix:** Re-run with `-D 500` (or `2000`) and read the DP of every site; a variant that already scores at itself (like GLA c.639+919G>A, 0.30) is found at `-D 50`. CI-SpliceAI is a second opinion only.

### SpliceAI: Tissue Agnosticism

**Trigger:** Variant in a tissue-specific gene (NEFM in neurons, MAPT brain, DMD muscle isoforms).

**Mechanism:** SpliceAI is trained on aggregate GENCODE annotation; tissue-specific events with weak constitutive use score low.

**Symptom:** Tissue-specific pathogenic variant has low SpliceAI delta; functional impact still observed in target tissue.

**Fix:** Use Pangolin per tissue or SpliceTransformer flags; require RNA validation in disease-relevant tissue.

### Pangolin: Masking and Out-of-Training Tissue

**Trigger:** `-m True`, or a disease tissue not represented in Pangolin's 4-species, 4-tissue (Cardoso-Moreira 2019 developmental) training set.

**Mechanism:** The mask erases real gains and losses (table above); Pangolin extrapolates poorly outside its training tissues.

**Symptom:** A splice-site loss or pseudoexon gain visible at `-m False` is 0.00 at `-m True`; score uncalibrated for the queried tissue and disagrees with patient RNA-seq.

**Fix:** Run `-m False`; fall back to SpliceAI for tissues not in the Pangolin training set, or run patient RNA-seq directly.

### MMSplice: Atypical Events

**Trigger:** Variant affecting a non-cassette event (MXE, complex multi-junction, AFE/ALE), or a deep-intronic variant.

**Mechanism:** MMSplice modular model is trained primarily on cassette exon events and scores variants only near an annotated exon.

**Symptom:** MMSplice ΔPSI doesn't match other predictors or empirical data, or the variant has no row.

**Fix:** Use SpliceAI for non-cassette events; restrict MMSplice to cassette exon contexts.

### CADD-Splice: Loss of Component Information

**Trigger:** Wanting to know which sub-component drove a high CADD-Splice score.

**Mechanism:** CADD-Splice combines splice predictors with the CADD annotations into a single C-score; sub-component contributions are abstracted.

**Symptom:** "High CADD-Splice score but unclear why."

**Fix:** Run SpliceAI and MMSplice separately to see which contributed.

### Branchpoint Variants: Low Information Motif

**Trigger:** Variant in the BPS region (-18 to -40 from 3'ss).

**Mechanism:** BPS motif (yUnAy) has low information content; CNNs struggle to learn the consensus.

**Symptom:** Confirmed BPS variant scores SpliceAI delta <0.2 despite functional disruption.

**Fix:** Require RNA validation (see Branchpoint Variant Detection).

## Population Database Lookup

| Database | Use for |
|----------|---------|
| gnomAD v4 | Allele frequency; SpliceAI annotations integrated |
| ClinVar | Existing classifications; SpliceAI integrated since 2020 |
| SpliceVarDB | Curated splice variants with experimental RNA validation |
| dbNSFP4 | Pre-computed splice scores aggregated |
| Recount3 | Tissue-specific PSI lookups from public RNA-seq |
| GTEx sQTL v8 | Tissue-specific splicing QTLs across 49 tissues |
| MaveDB | Splice MAVE results (e.g. BRCA1 saturation; Findlay 2018 *Nature*) |

Always check ClinVar first for existing classifications; cross-reference with gnomAD for population frequency before committing to PP3/PP4.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `ImportError: pkg_resources` / `No module named 'pkg_resources'` (SpliceAI, MMSplice, Pangolin) | setuptools >= 81 | `pip install "setuptools<81"` |
| SpliceAI/Pangolin skip a record: `Skipping record (ref issue)` / `Mismatch between FASTA (ref base: G) and variant file (ref base: A)`, exit 0 | REF does not match the FASTA (wrong build or coordinates) | Assert REF against the FASTA before scoring; report records missing from the output |
| Pangolin `Mismatch between FASTA (ref base: c) ...` on a correct REF | Soft-masked (lower-case) FASTA, e.g. UCSC `hg38.fa`; the comparison is case-sensitive | Use an upper-case FASTA (GENCODE genome) or `awk '/^>/{print;next}{print toupper($0)}'` |
| `ValueError: Duplicate ID CDS:ENST...` from `gffutils.create_db` | Real GENCODE GFF3 has repeated IDs | Use Pangolin's `create_db.py` on the GTF |
| `ModuleNotFoundError: pyfastx` / `vcf` / `torch` from Pangolin | `setup.py` declares no dependencies | Install them (see Install); `TypeError: Info.__new__() missing 'type_code'` -> `PyVCF3==1.0.0` |
| `pangolin: WARNING, skipping variant: Variant not contained in a gene body` | Contig names differ from the DB or gene absent from the (canonical) DB | Match chromosome naming; rebuild the DB from the full GTF |
| `Variant format not supported` (Pangolin) | REF and ALT both longer than one base with different lengths, or an allele with no ACGT (N-only) | Decompose into simple SNV/indel records; SpliceAI writes `.` scores for N-only alleles |
| SpliceAI `OSError: Can't write record` | Symbolic ALT such as `<DEL>` | Remove or expand structural variants |
| `mmsplice: Input sequence ... intron length cannot be longer than the input sequence` | FASTA truncated or not matching the GTF | Check FASTA size/index against the reference; use the GENCODE basic GTF |
| MMSplice output has no row for a variant | Variant not near an annotated exon (silent) | Compare output IDs with the input VCF |
| `SpliceVault` returns no rows | Variant not a catalogued splice-site SNV, or wrong REF/ALT/build | Check with SpliceAI which site is affected; GRCh38 only |
| `VariantValidator: invalid HGVS` | Wrong reference transcript or build | Specify NM_*.* version explicitly |
| A `chr`-prefixed VCF against a bare-contig FASTA (`chrX` vs `X`) | Naming mismatch | Both SpliceAI 1.3.1 and Pangolin 1.0.2 scored such records on the audit data; still keep names consistent (`bcftools annotate --rename-chrs`) if a tool skips |

## Common Pitfalls

- **Using SpliceAI score alone for clinical reporting** — combine with concordant predictors and ideally RNA validation; see ClinGen rules above.
- **Reading a missing score as benign** — `.`, no `SpliceAI=` tag, no Pangolin/MMSplice row = not scored.
- **Forgetting NMD direction** — confirmed splice disruption needs NMD-status check. Last-exon PTCs escape NMD and can be dominant-negative or gain-of-function.
- **Trusting LLMs for variant interpretation** — use as orchestrators on top of SpliceAI/VariantValidator/ClinVar; all clinical-grade calls require human expert sign-off.
- **Skipping HGVS validation** — invalid HGVS leads to silent reference-transcript mismatches; always run VariantValidator first.
