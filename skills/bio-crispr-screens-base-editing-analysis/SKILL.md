---
name: bio-crispr-screens-base-editing-analysis
description: Analyzes base-editing screens for variant function. Covers library design (Hanna 2021 ClinVar-scale CBE screen benchmarked on BRCA1/2, Cuella-Martin 2021 DDR saturation), CBE vs ABE chemistry choice (BE3/BE4 vs ABE7.10/ABE8.20/ABE8e), editing-window math (positions 4-8 from PAM-distal end; 4-7 for ABE7.10), bystander-edit quantification and the variant-call ambiguity it creates, sgRNA-efficiency filtering before hit calling, indel byproduct interpretation, the substitution-vs-indel diagnostic, variant annotation against ClinVar / COSMIC, and the Broad be-validation-pipeline. Use when designing a BE variant screen, choosing CBE vs ABE for a specific edit, interpreting bystander-confounded hits, distinguishing functional signal from indel artifact, integrating CRISPResso2 output with screen scoring, or deciding BE vs PE for SNV installation.
tool_type: mixed
primary_tool: CRISPResso2
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked against CRISPResso2 2.3.4 (2026-09-16, Docker `pinellolab/crispresso2:latest`) and BE-Hive git HEAD (maxwshen/be_predict_bystander, 2026-09-16) real output -- the parsers in this Skill assume that output schema, not the file layouts described in older CRISPResso2 docs. Also tested with pandas 2.2+, biopython 1.83+, numpy 1.26+, scipy 1.12+, scikit-learn 1.4+; Broad be-validation-pipeline notebooks (repo HEAD).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `CRISPResso --version` (install: `conda install -c bioconda crispresso2`)
- Python: `pip show CRISPResso2`; BE-Hive is a GitHub clone (maxwshen/be_predict_bystander), not a PyPI package

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Base Editing Screen Analysis

**"Analyze my base-editor variant-function screen"** -> Quantify per-sgRNA target-base conversion, bystander rate, and indel byproducts from amplicon sequencing; filter on editing efficiency; map each sgRNA to its intended SNV (target + bystander pattern); compute per-variant fitness from the screen log-fold change; reconcile target vs bystander variant attribution; annotate against ClinVar / COSMIC.

- CLI: `CRISPResso --base_editor_output` for per-amplicon BE quantification
- CLI: Broad `be-validation-pipeline` for end-to-end pooled-screen analysis with editing-efficiency filtering
- Python: `BE-Hive` (Arbab 2020) for editing-efficiency prediction; clone maxwshen/be_predict_bystander and import via sys.path -- see `references/be-hive-prediction.md`
- Web: `BE-Designer` (Hwang 2018, RGEN Tools) for variant-encoding sgRNA design

## Reference Files

| Need | Read |
| --- | --- |
| Predict per-spacer editing efficiency / bystander outcomes with BE-Hive (50nt substrate) | `references/be-hive-prediction.md` |
| Tile spacers across a protein region, annotate target vs bystander edits | `references/library-design.md` |
| After the screen: filter sgRNAs by editing efficiency (>50% convention), attribute signal to target vs bystander, aggregate sgRNA scores to variants | `references/screen-analysis.md` |
| Hanna 2021 / Cuella-Martin 2021 design and results | `references/published-screens.md` |
| Post-process CRISPResso2 BE amplicon output with the Broad notebooks | `references/be-validation-pipeline.md` |

Runnable code lives in `scripts/` (each has a header with inputs and a usage line; the reference files above show the invocation): `behive_predict.py`, `find_be_spacers.py`, `filter_by_editing_efficiency.py`, `deconvolute_bystander.py`, `aggregate_variant_scores.py`. `examples/base_editing_analysis.sh` has the CRISPResso2 CLI runs.

## Base Editor Chemistry Selection

| Editor | Reaction | Editing window | Indel byproduct rate | When to use |
|--------|----------|----------------|----------------------|-------------|
| BE3 (Komor 2016) | C->T (also G->A on opposite strand) | Pos 4-8 from PAM-distal end | 5-10% | Original; superseded |
| BE4 / BE4max (Koblan 2018) | C->T | Pos 4-8 | <5% | CBE standard |
| eA3A-BE3 | C->T narrow specificity | Pos 5-7 | <5% | Specifically TC contexts (eA3A prefers TC) |
| ABE7.10 (Gaudelli 2017) | A->G (T->C opposite strand) | Pos 4-7 | <2% | First ABE; slow at non-TA contexts |
| ABE8.20 (Gaudelli 2020) | A->G | Pos 4-8 | <2% | Modern ABE; high activity |
| ABE8e (Richter 2020) | A->G | Pos 4-8 | <2% | Highest editing activity; more processive than ABE7.10 |
| evoCDA-BE | C->T (broader) | Pos 1-9 | 5-10% | Larger editing window; more bystander |
| CGBE1 (Kurt 2021) | C->G | Pos 5-7 | 5-10% | C-to-G transversion; rare use |
| GBE (Zhao 2021) | C->G or C->A | Pos 4-7 | 5-10% | Transversions; less mature |

**Decision rule:** For a target SNV at position 4-8 of a candidate spacer with no bystander Cs/As in the same window, BE3-BE4 or ABE7.10 is sufficient. For high-throughput variant scanning where bystander tolerance must be minimized, use eA3A-BE3 (TC contexts only) for C->T, or ABE7.10 rather than ABE8e/ABE8.20 for A->G -- its 4-7 window is the narrowest ABE.

## Editing Window Math

**Why this matters for postdoc-level use:** Base editors are tethered to dCas9 (or nCas9) and the deaminase acts on the displaced ssDNA "R-loop" formed when Cas9 binds. The deaminase has a fixed reach -- positions 4-8 from the PAM-distal end of the protospacer for canonical BE3/BE4, and 4-7 for ABE7.10. Outside this window, editing efficiency drops by 10-50x.

```
PAM-distal end                                                            PAM-proximal
   |                                                                          |
   1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20    NGG
                  ^^^^^^^^^^^
                  Canonical editing window (positions 4-8)

   For BE4max: positions 4-8 are 5-50x more efficient than positions 1-3 or 9-13 (ABE7.10: 4-7)
   For SpABE8e: positions 4-8 (Richter 2020), matching the corresponding CBEs rather than ABE7.10's narrower 4-7
   For evoCDA-BE: window 1-9 (broader; more bystander)
```

**Critical implication for variant interpretation:** If the intended edit is at position 5 and there is an additional editable C/A at position 7, both will be edited in the same molecule. The screen scores the *combination* of edits, not the intended one alone. This is bystander confounding.

## Cas9 vs Base Editor vs Prime Editor for Variant Installation

| Approach | What it does | Bystander | Indels | When to use |
|----------|--------------|-----------|--------|-------------|
| Cas9 + HDR template | Installs precise edit + template | None | High (NHEJ competition) | When precise edit needed; high indel byproduct |
| Cas9 (no template) | Random indels at cut site | None | 70%+ | Loss-of-function; not variant-specific |
| CBE (BE3/BE4) | C->T at editing window | Yes (multiple Cs) | <5% | C->T variants with manageable bystanders |
| ABE (ABE7.10/ABE8e) | A->G at editing window | Yes (multiple As) | <2% | A->G variants; clean for single-A spacers |
| CGBE / GBE | C->G or C->A | Yes | 5-10% | Transversions; rare use cases |
| Prime editor (PE2/PE3) | Templated edit; any base change | None | 1-3% | Precise variants; lower efficiency |

**Decision:** For C->T or A->G with available editing window: base editor is preferred (higher efficiency than PE). For other transitions/transversions, multi-base edits, or zero-bystander requirements: prime editor.

## Validation Strategy

| Tier | Validation requirement |
|------|-------------------------|
| Tier 1 (high confidence) | BE + PE concordant at same variant + arrayed confirmation (convergent BE + PE is the gold standard for pathogenicity calls) |
| Tier 2 (medium) | BE alone, multiple sgRNAs converge despite bystander differences |
| Tier 3 (exploratory) | Single sgRNA hit; bystander confounded; not interpretable |

## Failure Modes

### Mostly indels in BE sample

**Trigger:** Cas9 contamination, wrong vector (e.g., used pCas9-BE3 plasmid but selected on Cas9 line), or evoCDA-BE / broader-window chemistry.
**Mechanism:** Cas9 cuts dsDNA; BE relies on nicked-ssDNA deamination. Cas9 expression in the same cell creates indels.
**Symptom:** Substitution-vs-indel ratio <3 in CRISPResso output.
**Fix:** Verify vector (nCas9-BE3 not Cas9-BE3); confirm cell line lacks Cas9 background; restrict to specifically engineered BE-cell lines.

### High editing but no biological signal

**Trigger:** Bystander C/A is dominating; intended variant is not the perturbation driving phenotype.
**Mechanism:** When target is at position 5 and bystander is at position 7, the molecule carries both; phenotype is from the bystander.
**Symptom:** Strong screen signal but variant attribution unclear.
**Fix:** Run orthogonal prime-editor scan of the same intended variants; restrict library to bystander-free spacers when possible; deconvolute via allele-frequency table.

### sgRNA shows perfect editing but no fitness signal

**Trigger:** Intended variant is silent or compensatory; the protein function is unchanged.
**Mechanism:** Variants can be tolerated; not all variants are LoF or GoF.
**Symptom:** High editing efficiency (>70%) but per-sgRNA LFC near zero.
**Fix:** Expected outcome for many variants; flag silent / compensatory variants in the report.

### Low editing across all guides

**Trigger:** Wrong cell line for the BE; cell line has poor BE activity (some lines lack APOBEC or have low expression).
**Mechanism:** BE efficiency depends on cell-line expression of TadA or APOBEC components.
**Symptom:** Median editing <30% across library.
**Fix:** Test in a BE-validated cell line (HEK293T, U2OS, K562 generally work); pilot before full screen.

### Library missing intended-variant sgRNAs

**Trigger:** No NGG-adjacent spacer places target base in editing window for that codon.
**Mechanism:** Editor window is fixed; some codons cannot be targeted with given chemistry.
**Symptom:** Specific variants absent from screen.
**Fix:** Use PAM-relaxed BE variants (SpRY-CBE, SpRY-ABE); use prime editor for variants outside BE accessibility; accept that some variants cannot be installed.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Editing window | Positions 4-8 from PAM-distal end (BE3/BE4); 4-7 (ABE7.10); 4-8 (SpABE8e) | Komor 2016; Gaudelli 2017; Richter 2020 |
| Editing efficiency for screen power | >30% (primary); >50% (validation) | Field convention (BE variant screens) |
| Indel byproduct (clean BE) | <5%; <2% for ABE | Koblan 2018 (BE4max); Gaudelli 2017 (ABE) |
| Substitution-vs-indel ratio | >10 (clean BE); <3 (Cas9-like) | CRISPResso2 diagnostic |
| Bystander rate (target attribution) | <10% acceptable; <5% ideal for clean attribution | Application-dependent |
| Cell-line BE activity (pilot) | >30% editing at validated target | Below = wrong cell line for BE |
| Per-amino-acid sgRNA density | 10-15 (saturation designs); 5-8 (smaller screens) | Tradeoff with library size |
| Plasmid pool evenness | Gini <0.1 | Verify by sequencing the pool before packaging |
| Lentiviral MOI | 0.3 | One integrant per cell |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Substitution-vs-indel ratio <3 | Cas9 contamination or wrong BE | Verify vector / cell line; pilot first |
| All edits at bystander positions | Target base outside window | Re-design spacer with target at pos 4-8 |
| Variant attribution unclear | Bystander confounding | Run orthogonal PE; restrict library |
| Library lacks intended variant | No NGG-PAM accessibility | SpRY-CBE; prime editor; accept exclusion |
| Editing <30% library-wide | Cell-line BE inactivity | Re-validate cell line |
| Hit list dominated by single sgRNA | Bystander-driven phenotype | Cross-check with bystander-free sgRNAs |

## References

- Komor AC et al. 2016. *Nature* 533:420. BE3.
- Gaudelli NM et al. 2017. *Nature* 551:464. ABE7.10.
- Koblan LW et al. 2018. *Nat Biotechnol* 36:843. BE4max + improved CBE.
- Richter MF et al. 2020. *Nat Biotechnol* 38:883. ABE8e; phage-assisted evolution of ABE7.10.
- Lapinaite A et al. 2020. *Science* 369:566. ABE8e mechanism.
- Gaudelli NM et al. 2020. *Nat Biotechnol* 38:892. ABE8 series (ABE8.20).
- Hanna RE et al. 2021. *Cell* 184:1064. Massively parallel BRCA1/2 variant function via CBE.
- Cuella-Martin R et al. 2021. *Cell* 184:1081-1097. CBE saturation across 86 DDR genes (BRCA1/2 plus others).
- Arbab M et al. 2020. *Cell* 182:463. BE-Hive prediction of editing outcomes.
- Anzalone AV et al. 2019. *Nature* 576:149. Prime editing (PE2/PE3).
- Clement K et al. 2019. *Nat Biotechnol* 37:224. CRISPResso2.
- Kurt IC et al. 2021. *Nat Biotechnol* 39:41. CGBE1.

## Related Skills

- crispr-screens/crispresso-editing - CRISPResso2 BE/PE mode and allele tables
- crispr-screens/library-design - base-editor library design
- crispr-screens/prime-editing-screens - Orthogonal PE for variant attribution
- crispr-screens/hit-calling - Variant-level hit aggregation
- crispr-screens/screen-qc - Editing-efficiency QC
- crispr-screens/drugz-chemogenomic - drugZ for BE drug-modifier screens
- clinical-databases/clinvar-lookup - Variant pathogenicity annotation
- variant-calling/variant-annotation - VEP for predicted amino acid changes
