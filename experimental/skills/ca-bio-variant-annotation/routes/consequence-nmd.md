# Consequence and NMD (inputs to PVS1)

No script. Use Sequence Ontology terms (`stop_gained`, `frameshift_variant`, `splice_donor_variant`, `splice_acceptor_variant`, `start_lost`, `stop_lost`, `missense_variant`, `inframe_deletion`, `synonymous_variant`), which VEP emits natively. SnpEff and ANNOVAR map to mostly equivalent terms but differ on splice-region width and finer intronic terms, so term-level matching across engines is impossible.

Whether a null consequence supports PVS1 depends on facts the impact bucket does not know:

- Premature stop more than about 50-55 nt upstream of the last exon-exon junction triggers nonsense-mediated decay: strong loss of function.
- A stop in the last exon, or within about 50 nt of the final junction (3' end of the penultimate exon), escapes NMD; the truncated protein is made. Downgrade PVS1 (Strong, Moderate or Supporting) by how much protein or which domains are lost.
- Single-exon genes have no junctions, so a nonsense variant escapes NMD by construction. A stop near the start can reinitiate downstream.
- PVS1 also needs loss of function to be the established disease mechanism. A null allele in a gain-of-function or dominant-negative gene must not fire it.
- Splice variants: the predicted result (exon skip, intron retention, cryptic site) determines PVS1, not the donor/acceptor term alone; SpliceAI does not report the result (`routes/predictors.md`).
- "HIGH impact" is never "PVS1 met". Combining into a classification is `variant-calling/clinical-interpretation`.

Done when each null variant has its exon position, NMD status and the gene mechanism stated.
