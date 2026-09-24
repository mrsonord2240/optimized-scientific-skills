# Motif Analysis with the Correct Background

**Goal:** Find kinase/writer motifs around the modified residue without rediscovering amino-acid composition bias.

**Approach:** Use the `Sequence window` (+/-15 residues, 31-mer) MaxQuant already provides, centered on the site. The background MUST be an experiment-matched S/T/Y set drawn from the identified proteins (or a central-residue-preserving shuffle), NOT the whole proteome or IUPAC-random -- those just report the composition of phospho-rich disordered regions. motif-x and MoMo p-values are only valid when the background is built this way.

```bash
python scripts/motif_enrichment.py --sites 'Phospho (STY)Sites.txt' --fasta search.fasta --out motif_enrichment.csv
# foreground = the regulated-up class-I sites only: add --foreground-ids up_sites.txt (one Gene_S473 id per line); --residues Y for tyrosine
```

`scripts/motif_enrichment.py` also exposes `matched_background()` and `motif_enrichment()` for import. It builds the foreground from the `Sequence window` of the class-I sites (Reverse/contaminant dropped, `Localization prob` >= 0.75), the background from every S/T window of the phospho-identified proteins in the FASTA (central residue preserved), then a one-sided Fisher test per (offset, residue) with BH across all tests.

For a publication-grade enrichment logo, hand the foreground and a matched background to a dedicated tool (motif-x / MoMo) and render with data-visualization/sequence-logos.
