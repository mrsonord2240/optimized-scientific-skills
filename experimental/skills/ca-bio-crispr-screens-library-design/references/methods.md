## Library Chemistry Decision Tree

| Goal | Chemistry | Canonical library | Guides/gene | TSS / target window |
|------|-----------|-------------------|-------------|---------------------|
| Loss-of-function essentiality, fitness | SpCas9 KO | Brunello, TKOv3, Avana | 4 (Brunello), 4 (TKOv3), 6 (Avana) | Constitutive exons, prefer aa 5-65% from N-terminus |
| Knockdown of non-cuttable genes, dosage-sensitive | dCas9-KRAB (CRISPRi) | Dolcetto, Horlbeck v2 | 6 (Dolcetto), 5 (Horlbeck) | Optimum +25 to +75 downstream of FANTOM5 TSS, searched out to -50/+300 (Dolcetto, Sanson 2018); -25 to +500 (Horlbeck v2) |
| Gain-of-function, gene activation | dCas9-VP64 / SAM / SunTag (CRISPRa) | Calabrese, Horlbeck-CRISPRa | 6 (Calabrese), 5 (Horlbeck) | -150 to -75 from TSS (Calabrese); -550 to -25 (Horlbeck v2) |
| Paralog buffering, GI screens | enAsCas12a multiplex | Inzolia, in4mer | 4-guide arrays | Constitutive exons |
| Variant function, SNV scanning | CBE / ABE | Custom tiling library | Tile editing windows | Editing window pos 4-8 from PAM-distal end |
| Precise edit, indel-free | Prime editor | Custom PRIDICT-designed | Tile pegRNAs | Anywhere with NGG PAM within 30 nt of edit |

Window functions for the CRISPRi/a rows: `references/crispri-crispra-tss.md`. Library sizes and PAM/enzyme choices: `references/library-catalog-and-pam.md`.

**Off-the-shelf first:** for a new screen in a well-characterized cancer line, use the published Brunello / Dolcetto / Calabrese / Inzolia pool (Addgene) rather than re-designing it; you inherit the community validation and the calibration of MAGeCK / BAGEL2 / Chronos against that library. Cas9 and Cas12a libraries are not interchangeable (different enzyme and PAM); do not mix them in one screen.

**Fails when:**
- CRISPRi/a targeting wrong TSS: any TSS without FANTOM5 CAGE evidence is suspect; guides positioned against the wrong TSS lose most of their knockdown.
- Cas9 KO of essential paralogs: single-KO buffering hides paralog-redundant essentials (42% of constitutively expressed genes never score, Dede 2020); switch to Cas12a multiplex.
- Base editor over an exon-intron boundary: editing-window bystanders create splice variants instead of the intended SNV.

## On-Target Scoring: Algorithmic Taxonomy

| Predictor | Year | Training set | Strengths | Fails when |
|-----------|------|--------------|-----------|------------|
| Doench Rule Set 1 | 2014 | Flow-sorted GFP+ knockouts | Simple, interpretable | Limited training data; sub-optimal at >NGG context |
| Doench Rule Set 2 / Azimuth 2.0 | 2016 | 1,841 flow-cytometry guides (Doench 2014) plus new guides tiling additional genes | Gold-standard for SpCas9; basis of Brunello | Trained on dropouts; under-predicts efficacy for nuclear-localized targets |
| DeepSpCas9 | 2019 | 12,832 synthetic targets integrated in HEK293T | Spearman ~0.77 vs measured indel frequency on held-out data | Black-box; sensitive to chromatin/context features it wasn't trained on |
| CRISPRon | 2021 | High-throughput indel sequencing | Best for therapeutic-grade target nomination | Slow per-guide; over-fits to its specific cell line |
| DeepHF | 2019 | ~171k guides in HEK293T (WT 55,604; eSpCas9(1.1) 58,167; HF1 56,888) | Separate model per enzyme variant, including WT | Pick the model matching the enzyme actually used |

**Reconciliation:** When predictors disagree, prefer the model whose training cell line matches the screen line (DeepSpCas9 was trained on synthetic targets integrated in HEK293T). For Brunello selection, Azimuth/Rule Set 2 is sufficient because the library was built with it -- introducing a different scorer creates apples-to-oranges ranking with the original library.

## Off-Target Scoring

| Score | Year | Math | Cutoff convention |
|-------|------|------|--------------------|
| MIT (Hsu) | 2013 | Position-weighted mismatch penalty | Specificity score 0-100, higher is better; CRISPOR treats >=50 as a good guide |
| CFD (Doench) | 2016 | Position+nucleotide-specific penalty fit on Brunello | Per-site CFD >0.2 counts a candidate off-target (Doench 2016); CRISPOR's aggregate CFD specificity score is 0-100, higher is better |
| Elevation | 2018 | ML on CFD + mismatch positions | Tighter than CFD |

CFD remains the default for genome-wide library design. **Critical pitfall:** CFD penalizes only mismatches, not bulges; for ≤1 mismatch + 1-bp bulge off-targets, validate empirically with GUIDE-seq or CIRCLE-seq. CRISPOR reports both the MIT (Hsu) and CFD guide specificity scores in a single output.

## Score and Rank sgRNAs for a Target Gene

**Goal:** Generate ranked sgRNA candidates for a single gene, jointly scored on on-target activity (Rule Set 2 / Azimuth) and off-target liability (CFD).

**Approach:** Identify all PAM-adjacent 20-nt protospacers in the target gene's coding sequence, retain only those in the first 5-65% of the protein (constitutive-exon convention from Brunello), filter on GC 30-70% and absence of poly-T (≥4 Ts terminates U6), score on-target (default: this Skill's GC heuristic; see Version Compatibility for the real-Rule-Set-2 alternatives) and off-target with CRISPOR, then **greedily select the top N guides that are also mutually independent** -- composition filters alone do not reject two candidates that overlap almost entirely (see `select_independent_guides` below).

The three functions live in `examples/design_library.py` (run `python examples/design_library.py`, or copy them):

- `find_sgrna_candidates(cds_sequence, pam='NGG', guide_length=20)` returns every PAM-adjacent protospacer on both strands with `spacer`, `strand`, `pos_in_cds`, `gc_frac`, already filtered for GC 30-70% and no `TTTT`. It uppercases the input first: the PAM/spacer regex matches uppercase ACGT only, so lowercase (soft-masked) input would silently return 0 candidates. The caller still filters by exon position and on-target/CFD score.
- `annotate_exon_position(candidates_df, cds_length)` keeps protospacers within the first 5-65% of the CDS (Brunello convention). N-terminal indels truncate the protein, very-N-terminal hits can be rescued by alternative initiation, and C-terminal hits miss functional domains (Doench 2016 Nat Biotech).
- `select_independent_guides(candidates_df, n_guides, min_spacing=5, score_col='score')` greedily picks the top-scoring candidates that are mutually >=`min_spacing` nt apart. The first two functions filter composition only, so nothing stops two candidates 1-4 nt apart (almost the same 20-nt spacer, the same cut site) from both counting toward the per-gene quota. Verified on a real gene (TP53, NM_000546.6): 12 candidates with no spacing filter included a pair 1 nt apart (19/20 nt shared); with the filter every pair is >=5 nt apart and the quota is still filled from the next-best candidates.

## Control Guides

A genome-wide library should include:

| Control type | Count | Purpose |
|--------------|-------|---------|
| Non-targeting (scrambled, no genomic match) | 500-1,000 (~1% of library) | Primary null distribution for CRISPRi/a; safe baseline for normalization |
| Safe-harbor (AAVS1, ROSA26-equivalent) | 50-100 | Cas9-only: absorbs cut-toxicity baseline (matters for amplicon-correction) |
| Olfactory receptors (presumed non-expressed) | 50-100 | Second null set for orthogonal normalization |
| Reference essentials (CEGv2 subset: e.g. RPS3, RPL11, EIF3A, POLR2A) | 50-100 | Internal positive control; QC dropout signal |
| Reference non-essentials (NEGv1 subset) | 50-100 | Internal negative control; BAGEL2 calibration |

**Critical pitfall:** Using only AAVS1 as the negative control in a Cas9 screen creates a normalization baseline biased toward "any cut is bad." Always add NTCs or non-essentials so that downstream median normalization and PR-AUC against CEGv2 work without baseline-shift artifacts. Too few NTCs (<100 in a 70k library) leave the null variance unstable: normalization and FDR rest on it, so gene-level p-values turn erratic and MAGeCK FDR fluctuates between runs.

## Library Composition for Specialized Screens

**Paralog buffering (Cas12a multiplex):** Build 4-guide arrays where positions 1-2 target gene A and positions 3-4 target paralog gene B. Inzolia covers ~4,435 paralog pairs within ~49k arrays. Singleton controls (gene A alone, gene B alone) must be included to score genetic interaction = double_KO_LFC - sum(single_KO_LFC).

**Base editor screens (tiling-library design):** Tile NGG-adjacent spacers across exons; ensure editing window (positions 4-8 from PAM-distal end) lands inside coding exons; flag bystander Cs/As in the window for downstream interpretation. Restrict to 50-90% editing efficiency a priori (filter out predicted low-efficacy guides) -- see bio-crispr-screens-base-editing-analysis.

**Tiling / regulatory dissection:** Dense (every 5-10 bp) CRISPRi or CRISPRa guides across the candidate region; CRISPRi has broader signal width (good for enhancer discovery) but Cas9-indel tiling has sharper resolution (good for pinpointing critical bases). Pair with CRISPR-SURF deconvolution.


## Library QC After Cloning

| Metric | Target | Failure mode if missed |
|--------|--------|------------------------|
| sgRNA detection (>25 reads/guide in plasmid pool) | ≥99% | Founder effect: missing guides cannot be screened; dropout impossible to distinguish from missing |
| Gini coefficient of plasmid pool | <0.1 | Synthesis defects or PCR bias; pool unfit for screening at standard 500x coverage |
| Skew ratio (top 10% / bottom 10%) | <2 (good), <5 (acceptable) | Skew >5 means underrepresented guides cannot generate statistical signal even at 1000x |
| % zero-count sgRNAs in plasmid pool | <0.5% | Plasmid bottleneck during cloning; re-amplify or re-clone |
| Replicate Pearson on plasmid pool (between sequencing technical replicates) | >0.99 | Sequencing artifact, not biology |

**Plasmid pool sequencing convention:** 200-500 reads per sgRNA before any biology (i.e. 15-40M reads for a 77k Brunello). This is the baseline against which all downstream depletion is computed; sequencing the plasmid is non-negotiable. If a metric misses its target, see `references/failure-modes.md`.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| GC content | 30-70% | Doench 2016 Nat Biotech: guides outside this range have low activity |
| Poly-T avoidance | ≤3 consecutive T | U6 Pol III terminator; ≥4 Ts terminates sgRNA transcription |
| Guides per gene (Cas9) | 4 (Brunello/TKOv3 standard); up to 6 (Avana, older) | Doench 2016 reports diminishing gene recovery below 4 sgRNAs/gene; returns flatten above 6 |
| Guide-to-guide minimum spacing | >=5 nt between any two selected guides for the same gene | Composition filters (GC/poly-T) alone don't reject near-duplicate protospacers; `select_independent_guides` enforces this so 4 "guides" are 4 independent cut sites, not fewer |
| MOI | 0.3 | Poisson: P(>=2 sgRNAs/cell) = 4% at MOI 0.3 |
| Coverage at infection | 500 cells/sgRNA | DepMap convention; 200x minimum, 1000x for noisy / in-vivo |

Windows: Library Chemistry Decision Tree. NTC counts: Control Guides. MIT/CFD cutoffs: Off-Target Scoring. Skew: Library QC After Cloning.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| sgRNA fails to express | Poly-T in spacer terminates U6 | Filter `TTTT` in design; this is the #1 silent failure |
| Hits include amplified loci (e.g. ERBB2 in HER2+) | Copy-number amplicon false-essentiality | See bio-crispr-screens-copy-number-correction |
| Paralog gene absent from hit list despite expression | Cas9 single-KO buffering | Switch to Cas12a multiplex; see bio-crispr-screens-combinatorial-screens |
| Cas12a oligo doesn't cut | Forgot Cas12a's TTTV PAM is 5' of spacer, not 3' | Re-orient: PAM-then-spacer for Cas12a, opposite of Cas9 |

No knockdown from a CRISPRi guide (wrong TSS): see `references/crispri-crispra-tss.md`. Plasmid-pool skew, synthesis dropouts, polyclonality: see `references/failure-modes.md`.
