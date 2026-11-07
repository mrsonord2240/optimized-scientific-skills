# i3_ora_background (not agent-visible)

## Pitfall and source
Over-representation analysis against a whole-genome background, when only a tissue-specific subset of genes is detectable, makes tissue-program gene sets look enriched in any DE list (Timmons et al. 2015 Genome Biol; Wijesooriya et al. 2022 PLoS Comput Biol). Design taken from that literature, not from any Skill.

## World (seed 20260919)
- 20,000 genes GENE00001-GENE20000 (`all_annotated_genes.txt`); 4,000 detected (`detected_genes.txt`).
- 60 sets in `pathways.gmt` (ids GS0001-GS0060, shuffled so id order carries no role):
  - 2 tissue-program sets (300 members each, 70% detected vs ~19% base rate). DE hits sit at chance level against the detected background (p = 0.31 each) but p = 7e-7, q = 1e-5 against the whole genome.
  - 2 true sets (90 and 120 members, 40% detected), 14 and 16 DE genes injected. q = 1.7e-7 and 3.3e-7 on the detected background (also significant on the others).
  - 56 nulls: 6 small (15-40 members; <10 detected, dropped by the size rule), 8 partly detected (38%), 5 under-detected (8%), 37 ordinary. No null has 8-12 detected members, so the 10-member filter has no borderline case.
  - ~25% of sets carry 1-3 `LOCnnnnnn` symbols absent from the annotation.
- DE list: 300 genes, all detected; 30 injected from the two true sets, remainder random detected genes. DE draw was rejection-sampled (attempt 209) so that: every null has p >= 0.02 on the detected background and on the genome background; tissue sets have p in [0.2, 0.95] (detected) and p <= 2e-5 (genome); >= 2 nulls have raw p in [0.02, 0.05).
- Truth = the two true set ids (GS0012, GS0048), equal to the output of the reference `correct()` (background = detected genes, hypergeometric, BH, size 10-500 after restriction).

## Margins on the correct background
Smallest q among non-true sets = 0.396 (smallest raw p 0.030); largest true q = 3.3e-7. Nothing is near 0.05.

## Solvers (python solvers.py)
- correct: passes.
- trap whole-genome: 4 sets (adds both tissue sets), fails.
- trap pathway-union background (8,029 genes, 152 DE in it): 4 sets, same adds, fails.
- trap raw p < 0.05: adds 2 nulls (GS0002, GS0050), fails.
- trap sets not restricted to the background (K from raw GMT size): 0 sets, fails.
- Variants that PASS: Fisher exact, no size filter, filter >10, filter on raw GMT size, test only sets with >=1 hit, Bonferroni.

## Tolerance rationale
Exact set match (Jaccard threshold 1). Effect sizes chosen so every variant agrees; verified above.

## Blind spots / ambiguity
- No blind spot found among named traps.
- The size rule is stated in the task ("after restricting the set to the genes in your background"), which mildly signals that a background choice exists, but does not say which. Both "detected" and "annotated" lists are plainly described, as a requester would.
- A solver that uses the detected list but does not restrict the sets to it (K = raw size) gets no hits; counted as a failure, arguably a plain mistake.
- A solver that uses genome background but drops sets by size after restriction gets the whole-genome answer; still fails.
- Pathway-union background is not spurious-free here because the union (8k genes) is 80%+ undetected genes; in a real GMT with mostly expressed genes it might behave closer to the correct background. That is a property of this fabricated world.
- Symbols are invented; nothing to recall. Real-world caveat: DE hits equal to detected genes by construction.
