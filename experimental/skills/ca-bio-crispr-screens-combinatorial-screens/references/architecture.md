## Combinatorial Architecture Decision Tree

| Goal | Architecture | Library | Why |
|------|--------------|---------|-----|
| Paralog buffering, identify synthetic lethal paralog pairs | enCas12a single-array 4-guide multiplex | Inzolia (Esmaeili Anvar 2024) | Cas9 single-KO misses paralog-buffered essentials (42% of constitutively expressed genes never score, Dede 2020) |
| Test specific pathway pair (e.g., DNA repair branches) | Big Papi (orthologous SaCas9 + SpCas9; two sgRNAs from U6 and H1 in pPapi) | Custom | Mature methodology; orthologous enzymes avoid repeated-element recombination |
| Combinatorial 3-way / 4-way knockout | in4mer (4-guide single Cas12a array) | Custom (in4mer) | Single transcript processed by Cas12a; multi-gene |
| Single-cell Perturb-seq with multi-pert per cell | Combinatorial Perturb-seq + Cas9 multiplex | Custom | Single-cell readout of multi-perturbation effects |
| Drug-modifier + KO interaction | Cas9 KO + drug treatment | Standard libraries | Drug as second "perturbation" |

**Fails when:**
- Dual-sgRNA constructs built from repeated U6/tracr elements: lentiviral recombination collapses them to a single perturbation
- Cas12a screens analyzed as Cas9 screens: MAGeCK normalization fails because Cas12a has different cut profile
- in4mer 4-guide arrays without all-singleton controls: GI scoring requires single-gene baselines

## Cas9 vs Cas12a for Multiplex

| Property | Cas9 paired (Big Papi) | Cas12a multiplex (in4mer / Inzolia) |
|----------|------------------------|--------------------------------------|
| Multiplex capacity per cassette | 2 sgRNAs (paired) | 4 (in4mer); 2 (standard Cas12a) |
| sgRNA processing | U6 and H1 promoters driving sgRNAs for two orthologous Cas9s | Single transcript processed by Cas12a itself |
| sgRNA inhibition with multiple targets | None | None (Cas12a's intrinsic processing handles all) |
| Library size for 1,000 pairs | ~4,000-6,000 paired cassettes (4-6 per pair) | ~2,000 arrays (2 per pair) plus singleton controls |
| Validated libraries | Limited (mostly custom) | Inzolia: ~49k 4-guide arrays covering 19,687 genes plus ~4,435 paralog pairs |
| Per-perturbation editing efficiency | High (each sgRNA independently) | Variable (Cas12a less efficient on some targets) |
| Best for | Pairwise GI of specific interest | Genome-scale paralog buffering; multi-gene perturbation |

**Recommendation:** For modern paralog screens, use Cas12a multiplex with the Inzolia library. It is ~30% smaller than a typical monogenic Cas9 library while additionally covering ~4,000 paralog pairs (Esmaeili Anvar 2024), which makes it more cost-effective at genome scale.

## The Paralog Buffering Phenomenon

**Dede et al 2020 *Genome Biol* 21:262** showed that a large share of constitutively expressed genes are never scored as essential in any Cas9 single-KO fitness screen (3,032 of 7,282; 42%), and that these never-essentials are strongly enriched for paralogs. The reason: gene paralogs perform redundant essential functions. Loss of one paralog is buffered by the other; only loss of both creates the essentiality phenotype.

**Quantified impact:** 24 synthetic-lethal paralog pairs identified in Dede 2020 across 3 cell lines; 19 of 24 (79%) reproduce in >=2 lines, 14 of 24 (58%) in all 3. These pairs were not findable by single-gene Cas9 screens, requiring combinatorial methodology.

**Examples:**
- **MAPK1 (ERK2) + MAPK3 (ERK1):** ERK family redundancy in proliferation
- **PIK3CA + PIK3CB:** PI3K alpha/beta redundancy
- **AKT1 + AKT2:** AKT family redundancy
- **HSP90AA1 + HSP90AB1:** HSP90 alpha/beta redundancy
- **STAG1 + STAG2:** Cohesion complex paralogs

Each is buffered: loss of one is tolerated; loss of both is lethal.


## Inzolia / in4mer 4-Guide Array Analysis

**Esmaeili Anvar 2024 *Nat Commun* 15:3577** introduced in4mer, a Cas12a multiplex architecture where each array contains 4 guides processed by Cas12a's intrinsic crRNA-processing activity. The Inzolia library is the canonical implementation, covering the protein-coding genome plus ~4,435 paralog pairs.

**Library design:**
- 4 guides per cassette (Cas12a single-transcript array)
- Per pair: 2 arrays carrying 2 guides per gene, with the guides presented in different order across the two arrays
- Includes singleton controls: each single gene is covered by 2 four-guide arrays (4 guides per gene, order swapped)
- ~49,000 total arrays covering 19,687 genes, ~4,435 paralog pairs, 376 triples, and 100 quads

```python
# Per-pair analysis from in4mer screen
def in4mer_pair_analysis(paired_counts_df, gene_pairs, value_cols):
    '''Aggregate cassette-level counts to per-pair statistics.
    paired_counts_df: rows = cassettes, with a cassette_id COLUMN (reset_index first if it is the index).
    gene_pairs: DataFrame with cassette_id and gene_A, gene_B columns.
    value_cols: the numeric sample/LFC columns to aggregate.
    '''
    merged = paired_counts_df.merge(gene_pairs, on='cassette_id')
    return merged.groupby(['gene_A', 'gene_B'])[value_cols].agg(['mean', 'std', 'count'])
```

## Failure Modes

### Dual-sgRNA construct recombines in the lentiviral vector

**Trigger:** A dual-sgRNA construct built from repeated elements -- two copies of the U6 promoter, or two copies of the SpCas9 tracrRNA scaffold.
**Mechanism:** Najm 2018 reports that repetitive elements in lentiviral vectors, including the U6 promoter and multiple copies of the tracrRNA sequence, drive high levels of recombination and reduce combinatorial screen efficiency. Big Papi avoids this by pairing two orthologous enzymes (SaCas9 + SpCas9), whose scaffolds differ, and expressing the two sgRNAs from distinct U6 and H1 promoters.
**Symptom:** Constructs collapse to a single perturbation; measured GI scores are diluted toward zero.
**Fix:** Use the pPapi architecture (orthologous Cas9s, U6 + H1) rather than duplicated U6/tracr elements; verify construct integrity by amplicon sequencing of clones.

### Cas12a screen with low editing efficiency

**Trigger:** Cas12a less efficient than Cas9 at some loci; some guides in the 4-guide array don't cut.
**Mechanism:** Cas12a editing rate varies by sequence context; some loci edit at <30%.
**Symptom:** Specific pairs missing expected effects despite cassette presence.
**Fix:** Pilot Cas12a efficiency at the loci before full screen; use enCas12a (enhanced) variant; for known low-efficiency loci, supplement with Cas9.

### GI scoring without singletons

**Trigger:** Library lacks single-gene controls (only paired knockouts).
**Mechanism:** GI = paired - expected_additive requires single-gene LFC; without them, expected cannot be computed.
**Symptom:** Cannot score GI; only paired LFCs available.
**Fix:** Design library to include singletons (place gene A with 3 placeholder guides; gene B with 3 placeholders); re-run with full design.

### Single-gene LFCs from different cell line

**Trigger:** Using public single-gene LFCs (e.g., DepMap) as the baseline for paired-screen GI scoring.
**Mechanism:** Single-gene effects are cell-line specific; using HCT116 single-gene LFCs to score K562 paired-screen GIs is invalid.
**Symptom:** GI scores look noisy; many false positives.
**Fix:** Include singleton controls in the screen; or use cell-line-matched DepMap data.

### Confounding cell-cycle / proliferation in GI scoring

**Trigger:** Paired KO of two cell-cycle-impacting genes; the double-effect saturates cell cycle.
**Mechanism:** If A_KO causes 50% growth arrest and B_KO causes 50%, the combined 75% arrest is already saturating proliferation; additive expectation overestimates double-effect, generating false "synthetic-rescue."
**Symptom:** GI scores positive for pairs of essential cell-cycle genes; biologically unexpected.
**Fix:** Use log-space (LFC) GI scoring rather than linear; saturation is less severe in log-space. Alternative: model with logistic / saturable response curve.

### Library skew amplifying noise

**Trigger:** Inzolia library has uneven cassette representation; some pairs at 10x lower coverage than others.
**Mechanism:** Standard library QC (Gini, skew) applies; low-coverage cassettes yield noisier LFCs.
**Symptom:** GI z-scores vary 2-3x across cassettes targeting the same pair.
**Fix:** Standard library QC; for low-coverage pairs, aggregate fewer cassettes but with more sequencing depth; or drop low-coverage pairs from analysis.

## Cross-Modality Validation

**Scope:** This Skill analyzes screen data for research purposes. A synthetic-lethal call here is a research lead, not a validated therapeutic target or a patient-treatment recommendation -- cell-line/screen-level synthetic lethality does not by itself establish patient-level clinical benefit. Any path from a screen hit to a therapy goes through the validation steps below and standard preclinical/clinical development; do not use screen output to make or imply a diagnostic or treatment decision for an individual patient.

For high-stakes synthetic-lethal hits (drug-target nomination), validate by:

1. **Orthogonal chemistry:** Re-validate with Cas9 if Cas12a, or vice versa
2. **Arrayed validation:** Single-knock-out arrayed setup with same cell line; quantify proliferation
3. **CRISPRi orthogonal:** Use dCas9-KRAB to confirm knockdown phenotype (no DNA damage)
4. **Pharmacological:** Inhibit paralog with drug; confirms target accessibility for drug development

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Synthetic lethal GI z-score | <-2 | Standard convention |
| Synthetic rescue GI z-score | >2 | Standard convention |
| No interaction | -1 to +1 | Within additive expectation |
| Cas9 paired-screen cassette count per pair | 4-6 | Standard library convention |
| Cas12a 4-guide arrays per paralog pair (Inzolia) | 2 (2 guides per gene, order swapped) | Esmaeili Anvar 2024 |
| Singletons in combinatorial library | At least 4-6 per single gene | For stable expected_additive |
| Cells per cassette for stable GI | 500+ at infection | Standard pooled-screen coverage |
| Cas12a editing efficiency for inclusion | >50% | Below = unreliable signal |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Dual-sgRNA construct acts as single | Recombination between repeated U6/tracr elements | Use pPapi (orthologous SaCas9 + SpCas9, U6 + H1) |
| Cas12a low editing | Locus-specific inefficiency | Pilot loci first; use enCas12a |
| Cannot compute GI | No singletons in library | Re-design to include all-singletons |
| GI scores noisy | Library skew | Standard library QC; aggregate cassettes |
| Many false "rescue" GIs | Saturation in linear-space | Use log-space (LFC) GI scoring |
| Too many SL/rescue calls at genome scale | Raw z-cutoff has no FDR control (see Multiple testing at genome scale, above) | Apply BH-FDR across all tested pairs' z-scores before calling hits |
| Drug-target paralog shows no GI in screen | Cell-line-specific buffering | Cross-validate with multiple lines |
