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
