---
name: bio-crispr-screens-prime-editing-screens
category: Data Analysis
description: Designs and analyzes pooled prime-editor (PE) screens for installing precise genetic variants without bystander confounding. Covers pegRNA design with PRIDICT and PRIDICT2 for predicting per-pegRNA editing efficiency, pegRNA architecture (spacer + scaffold + PBS + RTT), PE2/PE3/PE3b/PEmax variants, MOSAIC in situ saturation mutagenesis, the PRIME pooled-screen methodology (Ren 2023; ~3,699 ClinVar variant screens), chromatin context as a major locus-level determinant of PE efficiency, scaffold-incorporation and indel byproduct quantification with CRISPResso2, and the cross-modal validation strategy of PE + base-editor screens for variant function. Use when designing a pegRNA library for variant installation, choosing between BE and PE for a specific edit, predicting pegRNA efficiency before library synthesis, analyzing PE screen output, distinguishing intended-edit from scaffold-incorporation, or scaling PE screens to thousands of variants.
tool_type: mixed
primary_tool: PRIDICT2
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: PRIDICT2 git HEAD 2026-09-16 (https://github.com/uzh-dqbm-cmi/PRIDICT2), CRISPResso2 2.3.4 (via Docker, `pinellolab/crispresso2:latest`), ePRIDICT git HEAD ddcdba360c969469a536343dd89e29179961e367 (https://github.com/Schwank-Lab/epridict) with its `light` model, pandas 2.2+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `python pridict2_pegRNA_design.py single --help`; `python pridict2_pegRNA_design.py batch --help`
- CLI: `python epridict_prediction.py single --help` (run from the epridict repo root)
- Web: PRIDICT2 and ePRIDICT web interfaces at https://pridict.it/

ePRIDICT is Linux/macOS only -- `pybigwig` has no Windows wheel (use WSL on Windows).

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Prime-Editing Screen Analysis

**"Design or analyze a pooled prime-editor screen"** -> Design pegRNAs (spacer + scaffold + PBS + RTT) for intended edits, predict efficiency with PRIDICT2, filter pre-synthesis to efficient candidates, install variants in the screen, quantify intended-edit vs scaffold-incorporation vs indel via CRISPResso2, and aggregate to per-variant fitness scores.

- Python: `PRIDICT2` for pegRNA efficiency prediction
- Python: `ePRIDICT` for chromatin-context prediction; pair with PRIDICT2 rather than replacing it
- CLI: `CRISPResso --prime_editing_pegRNA_*` for amplicon-level analysis
- Workflow: pegRNA library design -> PRIDICT2 filtering -> ePRIDICT chromatin check -> screen execution -> CRISPResso2 quantification -> per-variant scoring
- Details: PRIDICT2 CLI and library filtering in `references/pridict2-batch-cli.md`; ePRIDICT install and prediction in `references/epridict-chromatin.md`; PRIME, MOSAIC and BE cross-validation in `references/pooled-screen-methods.md` (see Reference Files)

## Prime Editor Chemistry Comparison

| Editor | Year | Mechanism | Indel rate | Use when |
|--------|------|-----------|------------|----------|
| PE2 (Anzalone 2019) | 2019 | nCas9-RT fusion + pegRNA | 1-3% | Standard PE; lowest indel rate |
| PE3 | 2019 | PE2 + nick of opposite strand by additional sgRNA | 2-5% | Higher editing efficiency, slightly more indels |
| PE3b | 2019 | PE3 with edit-blocking ssgRNA | 1-3% | When PE3's added nick risks unwanted indels |
| PEmax (Chen 2021) | 2021 | Engineered RT + nCas9 | 1-2% | Higher editing rate per pegRNA |
| PE5max (Chen 2021) | 2021 | PE3 plus MMR inhibition (MLH1dn) on the PEmax architecture | 1% | Highest efficiency at favorable sites |
| PE6 / dual-pegRNA (2023) | 2023 | Engineered compact PE; twin-pegRNA systems | Variable | Specific applications |

**Decision rule:** For pooled screens at scale, PE2 or PEmax (single-guide architecture) is preferred over PE3, whose additional nicking sgRNA complicates library architecture. For specific high-efficiency edits, PEmax + PRIDICT2-optimized pegRNA.

## pegRNA Architecture

A pegRNA contains four critical elements that determine efficiency. **The 3' extension is
RTT-then-PBS, not PBS-then-RTT** -- PBS is the pegRNA's 3'-terminal element. Getting this order
backwards is not cosmetic: it makes real CRISPResso2 silently report 0% editing on a true 40%
sample (see the CRISPResso2 section below), because `--prime_editing_pegRNA_extension_seq` never
matches the edited allele under the wrong order.

```
5'  SPACER (20 nt)  -- standard sgRNA spacer; defines target locus via NGG PAM
    +
    SCAFFOLD (~80 nt) -- canonical or recoded scaffold (Chen 2021 recodes it to cut scaffold-incorporation byproducts)
    +
    RTT (Reverse Transcription Template, 10-30 nt) -- encodes intended edit; copied by RT.
        Spans the genomic region DOWNSTREAM of the nick (through the PAM and beyond) --
        this is the window PE actually replaces.
    +
    PBS (Primer Binding Site, 8-15 nt) -- the pegRNA's 3'-terminal element; primes reverse
        transcription off the nicked strand. Complementary to the protospacer sequence
        UPSTREAM of the nick -- it must never include PAM bases.
3'
```

Both RTT and PBS are written on the pegRNA as the **reverse complement** of their genomic
template regions (verified against the real PRIDICT2 CLI's own `RTrevcomp`/`PBSrevcomp`/`pegRNA`
columns and against real CRISPResso2 output -- see below).

**Worked example** (hand-derived, independently confirmed by a real CRISPResso2 2.3.4 run
recovering the planted 40.0% Prime-edited / 33.33% Modified-within-Reference exactly): spacer
`ACGTTGACCTGGAACGTTCA`, NGG PAM immediately 3' of it, C>T edit 11 nt downstream of the nick.
RTT genomic window (18 nt, downstream of nick, edit installed) = `TCATGGCGATCTGTAAGC`; its
reverse complement = `GCTTACAGATCGCCATGA`. PBS genomic window (13 nt, upstream of nick) =
`TGACCTGGAACGT`; its reverse complement = `ACGTTCCAGGTCA`. Extension = RTT-revcomp + PBS-revcomp
= `GCTTACAGATCGCCATGAACGTTCCAGGTCA` -- this is the value that belongs in
`--prime_editing_pegRNA_extension_seq`.

**Key design parameters:**
- **PBS length:** 11-13 nt typical; longer for high-GC contexts; PBS GC fraction critical (35-65% target)
- **RTT length:** 10-20 nt typical; longer for distant edits (10+ bp away from cut)
- **RTT-edit position:** intended edit at position 4-30 from cut site
- **Scaffold:** standard sgRNA scaffold OR the Chen 2021 recoded scaffold, which removes homology with the genomic target and cuts scaffold-derived byproducts. The recoding does not itself raise editing efficiency; that comes from MLH1dn and the PEmax architecture.

## CRISPResso2 for PE Quantification

```bash
# --quantification_window_size 25 widens the window to cover the edit
CRISPResso \
    --fastq_r1 pe_sample.fq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_spacer> \
    --prime_editing_pegRNA_spacer_seq <spacer> \
    --prime_editing_pegRNA_extension_seq <RTT-revcomp + PBS-revcomp, RTT first> \
    --prime_editing_pegRNA_scaffold_seq <scaffold> \
    --quantification_window_size 25 \
    --output_folder pe_results \
    --name sample_id

# Output: CRISPResso2 nests one level below --output_folder, as
#   pe_results/CRISPResso_on_sample_id/CRISPResso_quantification_of_editing_frequency.txt
# (i.e. <output_folder>/CRISPResso_on_<name>/; the commented filename alone is not the path).
# Prime-editing outcomes appear as extra amplicon ROWS (Reference / Prime-edited /
# Scaffold-incorporated), each with Unmodified%, Modified% and read counts.
# A correct extension puts the edited reads in the Prime-edited row; building it PBS-then-RTT
# instead yields Prime-edited NA/0 with those same reads counted as Reference-row substitutions
# (verified: 240/500 edited reads -> Prime-edited row with RTT-then-PBS, 0 -> with PBS-then-RTT).
```

## Failure Modes

### Low pegRNA efficiency despite high PRIDICT prediction

**Trigger:** Sequence-only prediction missed chromatin context.
**Mechanism:** Closed chromatin reduces Cas9 binding and RT activity; PRIDICT2 only sees sequence.
**Symptom:** PRIDICT2 predicts 60% efficiency; observed is 5%.
**Fix:** Score the target's chromatin context directly with ePRIDICT (K562 model; see
`references/epridict-chromatin.md`) rather than proxying it with a separate ATAC-seq track -- ePRIDICT
is the published companion to PRIDICT2 for exactly this gap. Flag pegRNAs whose ePRIDICT percentile is
low despite a high PRIDICT2 score; pilot those before committing library budget.

### High scaffold incorporation

**Trigger:** RTT too short relative to PBS, or RT processivity issue.
**Mechanism:** RT reads past edit into scaffold; resulting product is detectable but undesired.
**Symptom:** Scaffold incorporation >5%; intended edit efficiency low.
**Fix:** Re-design pegRNA with longer RTT; verify with PRIDICT2 score for scaffold_incorp; pilot at representative loci.

### PE2 cell line lacks RT expression

**Trigger:** PE2 construct expressed at low level; insufficient RT for productive editing.
**Mechanism:** PE2 requires high RT expression; some cell lines down-regulate.
**Symptom:** Library-wide editing <10%; not locus-specific.
**Fix:** Verify PE2 expression by Western blot; consider PEmax (higher activity); use better-validated cell lines (K562, HEK293T, U2OS).

### Multi-base intended edit but only one base installed

**Trigger:** Long RTT designed for multi-base edit; RT prematurely terminates.
**Mechanism:** RT processivity drops with longer RTT; multi-base edits often incomplete.
**Symptom:** Allele table shows partial-edit alleles (some bases installed, not all).
**Fix:** Re-design with shorter RTT covering only the closest edits; or use PE3 to nick opposite strand and force longer RT processivity.

### Library missing intended variant

**Trigger:** No suitable PAM/PBS/RTT combination for the intended edit.
**Mechanism:** PE requires NGG PAM within 30 nt of edit; rare edits cannot be installed.
**Symptom:** Specific variants absent from library.
**Fix:** Use SpRY-PE for relaxed PAM; accept that some variants cannot be PE-installed; consider BE if applicable.

### Wrong pegRNA-extension element order silently zeroes out real editing

**Trigger:** Building `--prime_editing_pegRNA_extension_seq` as PBS-then-RTT instead of RTT-then-PBS.
**Mechanism:** CRISPResso2 can never match the edited allele against the wrong-order extension sequence.
**Symptom:** Exit code 0; a true 40% Prime-edited sample reports 0.0% Prime-edited, with only a generic
"disproportionate percentages" / "substitutions outside quantification window" warning -- easy to miss.
**Fix:** Build the extension as RTT-revcomp + PBS-revcomp (see pegRNA Architecture above); if Prime-edited%
comes back near-zero on a library-wide basis, check element order before assuming a biological failure.

### Batch run exits 0 with an empty summary file

**Trigger:** Any batch run that produces zero successful pegRNA designs, or a bare `--summarize`.
**Mechanism:** A bare `--summarize` crashes argparse outright. Otherwise `summarize_top_scoring()` writes an
empty DataFrame whenever the output directory holds no per-sequence prediction CSVs, so three different
causes yield the byte-identical summary file `""` and exit code 0 (checked on PRIDICT2 git HEAD 2026-09-16;
none hangs or crashes): (1) wrong CSV header (`sequence` instead of `editseq`, prints
`Missing "editseq" column`); (2) correct header but zero data rows (prints `Designing pegRNAs for 0 sequences`);
(3) every variant lacks an NGG PAM within the search window (prints `No PAM (NGG) found in proximity of edit!`,
nick at most 25 bases from the edit).
**Symptom:** A "completed successfully" run whose summary file is 4 bytes, or a crash on `--summarize`.
**Fix:** Always pass `--summarize <K562|HEK>` and a CSV headed `sequence_name,editseq` with at least one row.
Read the stdout messages above, or check for per-sequence `predictions/<sequence_name>_pegRNA_Pridict_full.csv`
files, to find which cause it is; an empty summary is not evidence that PRIDICT2 found nothing. Variants with
no designable PAM cannot be PE-installed (see "Library missing intended variant").

## Cas9 vs BE vs PE for Variant Installation

| Approach | Bystander | Indels | Coverage | When to use |
|----------|-----------|--------|----------|-------------|
| Cas9 + HDR | None | High | Variable (depends on template integration) | Precise edits at scale; high indel byproduct |
| Base editor | YES | Low (<5%) | Limited by editing window | C->T or A->G at editable position |
| Prime editor | NONE | Low (<3%) | NGG-PAM within 30 nt of edit | Precise variants; multi-base; transversions |
| Cas9 (no template) | NONE | 70%+ | Anywhere with NGG | LoF only; not variant-specific |

**Decision tree:**
- C->T or A->G at editing-window position: BE (higher efficiency than PE)
- Multi-base / transversion / out-of-window: PE
- LoF without specifying variant: Cas9
- Random insertions: HDR (lower throughput than PE)
- iPSC / primary-cell variant: PE (no bystander confounding)
- Cancer-line variant scanning or drug-resistance variant: PE or BE; cross-validate both (chemistry-dependent) (see `references/pooled-screen-methods.md`)

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| PRIDICT2 efficiency for library inclusion | >50% | Project-chosen cutoff; PRIDICT2 prescribes none |
| Intended edit % for screen power | >5%; >20% at favorable sites | Field convention |
| Scaffold incorporation | <2% (clean PE); <5% acceptable | Empirical |
| Indel byproduct | <3% (PE2); <5% (PE3) | Anzalone 2019; Chen 2021 |
| PBS GC content | 40-55% | PRIDICT2 |
| PBS length | 11-13 nt | PRIDICT2 |
| RTT length | 10-20 nt | PRIDICT2 |
| Edit position from cut | 1-30 nt | Anzalone 2019 |
| Cell line for PE | K562, HEK293T, U2OS validated | High RT expression |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low editing across library | Cell-line RT inactivity | Verify PE2 expression; switch to validated line |
| Scaffold incorporation >10% | RTT too short | Re-design with longer RTT |
| Partial multi-base edits | RT processivity limit | Shorter RTT or PE3 |
| PRIDICT predicts but observes much lower | Chromatin context | Score the locus with ePRIDICT (`references/epridict-chromatin.md`); flag low-percentile targets; pilot before library order |
| Library missing variants | No NGG PAM | SpRY-PE; BE alternative |
| CRISPResso2 reports ~0% Prime-edited on a library that should edit | pegRNA extension built PBS-then-RTT instead of RTT-then-PBS | Rebuild extension as RTT-revcomp + PBS-revcomp |
| PRIDICT2 batch: `--summarize: expected one argument` | Bare `--summarize` flag | Pass a value: `--summarize K562` (or `HEK`) |
| PRIDICT2 batch: summary file is `""` (exit 0) | Wrong header (`sequence` not `editseq`), zero data rows, or no NGG PAM near any edit | Read the run's stdout for the cause; see Failure Modes, "Batch run exits 0 with an empty summary file" |
| PRIDICT2 batch: `FileNotFoundError` referencing `input/<file>` | CSV not placed in the CLI's default `./input` directory | `mkdir input` and move the CSV there, or pass `--input-dir` explicitly |
| PRIDICT2 batch `--summarize`: `FileNotFoundError` on the output directory itself | `--output-dir` does not exist yet -- `--summarize` lists existing `.csv` files there before the run starts | `mkdir -p <output-dir>` before running with `--summarize` |
| PRIDICT2 batch `--summarize`: `ValueError: Output directory is not empty. Please move or delete existing .csv files` | `--output-dir` exists but already holds a `.csv` -- including a previous PRIDICT2 run's own per-sequence `_pegRNA_Pridict_full.csv` files, which the run itself writes there | Pre-creating the directory is necessary but not sufficient: a second `--summarize` run into the same `--output-dir` always fails. Point it at a fresh directory (or clear the `.csv` files first) |
| PRIDICT2 batch with `--cores` above 3 | The tool documents 3 as the maximum ("Maximum 3 cores to prevent memory issues") but passes `--cores` straight through unvalidated, so a larger number is accepted silently | Use `--cores 3` (also the default); do not raise it |
| PE concordant with BE on transitions, disagrees on transversions | PE handles transversions BE doesn't | Expected; trust PE |

## Reference Files

| File | Read when |
|------|-----------|
| `references/pridict2-batch-cli.md` | Running PRIDICT2 single or batch mode, loading its output columns, filtering and picking top pegRNAs per variant |
| `references/epridict-chromatin.md` | Installing and running ePRIDICT, or interpreting a locus whose chromatin context contradicts its PRIDICT2 score |
| `references/pooled-screen-methods.md` | Planning a PRIME-style pooled screen or MOSAIC saturation library, or intersecting PE hits with a parallel BE screen |

## References

- Anzalone AV et al. 2019. *Nature* 576:149. Original PE2/PE3 (foundational prime editing paper).
- Mathis N et al. 2023. *Nat Biotechnol* 41:1151. PRIDICT v1 deep-learning pegRNA prediction.
- Mathis N et al. 2025. *Nat Biotechnol* 43(5):712 (published online June 2024). PRIDICT2 + chromatin context (current state-of-the-art).
- Chen PJ et al. 2021. *Cell* 184:5635. PEmax + engineered RT.
- Hsu JY, Lam KC, Shih J, Pinello L, Joung JK 2024 bioRxiv (doi:10.1101/2024.04.25.591078). MOSAIC in situ saturation mutagenesis via prime editing.
- Ren X et al. 2023. *Mol Cell* 83:4633. PRIME pooled prime-editing screen (~3,699 ClinVar variants); variant-installation scale.

## Related Skills

- crispr-screens/library-design - pegRNA library design
- crispr-screens/base-editing-analysis - Orthogonal BE for variant attribution
- crispr-screens/crispresso-editing - CRISPResso2 PE mode and quantification
- crispr-screens/hit-calling - Per-variant hit calling
- crispr-screens/screen-qc - Editing-efficiency QC
- variant-calling/variant-annotation - Annotate edited variants
- clinical-databases/clinvar-lookup - Variant pathogenicity
