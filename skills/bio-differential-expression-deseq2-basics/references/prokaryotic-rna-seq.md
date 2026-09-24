# Prokaryotic RNA-seq with DESeq2

## Prokaryotic RNA-seq

- Non-spliced aligners (BWA-MEM, Bowtie2) -- no introns.
- Polycistronic operons cause read-through between adjacent genes; confirm gene boundaries in the GFF.
- rRNA depletion essential (80-95% rRNA without poly-A selection, which prokaryotes lack anyway).
- Median-of-ratios fails under stress (see failure mode above) -- use `controlGenes` or spike-ins.
- KEGG organism codes are strain-specific (e.g., `pae` for P. aeruginosa PAO1): `clusterProfiler::search_kegg_organism()`.
- Annotation comes from Prokka or Bakta GFF; Ensembl/biomaRt are eukaryote-only.
