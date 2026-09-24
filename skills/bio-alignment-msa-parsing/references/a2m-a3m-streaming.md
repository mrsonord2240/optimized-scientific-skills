# A2M / A3M and Streaming

Read when the input is A2M/A3M or too large to load whole.

## A2M / A3M Conventions

A2M (HMMER) and A3M (HHsuite, ColabFold) encode insert vs match columns via case (uppercase = match column residue, lowercase = insert). A3M does not pad inserts across sequences and must be reformatted to A2M before loading as a rectangular MSA. Match-only extraction: `examples/a2m_a3m_io.py`, which reads with `SeqIO.parse` because HMMER 3.4 `hmmalign --outformat a2m` writes UNPADDED rows (no `.` characters; `AlignIO.read(..., 'fasta')` raises "Sequences must all be the same length"), whereas HH-suite pads. To get a rectangle from HMMER A2M use `pyhmmer.easel.MSAFile(path, format='a2m')`. See `alignment/alignment-io` A2M / A3M Conventions section for the full character table, BioPython load pattern, and the `reformat.pl` reference-sequence pitfall.

## Streaming Large Alignments

For Pfam-scale streaming (multi-gigabyte Stockholm or A3M databases that exceed RAM), use `pyhmmer.easel.MSAFile` with `compute_weights(method='pb')` for in-flight Henikoff weighting. See `alignment/alignment-io` Streaming Large Stockholm Databases section for the full code pattern.
