# Library format conversion and OpenSWATH decoys

Moved from SKILL.md. Read when converting a library between DIA-NN, OpenSWATH and Spectronaut, or building an OpenSWATH TSV/TraML and generating decoys.

## Convert Library Formats

**Goal:** Move a library between DIA-NN, OpenSWATH, and Spectronaut conventions without silently corrupting RT, intensity, modification, or decoy content.

**Approach:** Conversion is renaming columns AND reconciling units, not a copy. Check RT units (iRT ~ -25..150 vs normalized 0-1 vs minutes), intensity scaling (relative vs absolute), and modification notation (UniMod:35 vs +15.9949 vs Oxidation). For OpenSWATH, generate decoys with OpenSwathDecoyGenerator -- a target-only library has no null.

Column mapping and iRT-unit assertion: `python scripts/spectronaut_to_diann.py --in spectronaut.tsv --out diann.tsv`
(renames `RelativeIntensity` to `LibraryIntensity`, `FragmentMz` to `ProductMz`, and so on; refuses input whose iRT is outside -50..200).

**OpenSwathDecoyGenerator's real input requirements (verified on OpenMS 3.5.0):**
`TargetedFileConverter` converts a TSV without complaint even when it is unusable. The transition
list needs ALL of: a literal `Annotation` column (e.g. `y3^1`), chemically real theoretical fragment
m/z (not placeholders), and a per-precursor grouping column -- `transition_group_id`, unique per
PeptideSequence + PrecursorCharge (`FullUniModPeptideName` also works). Without the grouping
column, distinct peptides sharing a charge collapse into one `<Peptide id="_2">` group with no error
(2 targets become 1), and the decoy step then fails or reports wrong counts. Missing `Annotation` or
placeholder m/z gives `Number of decoy peptides: 0`.

`python scripts/build_openswath_tsv.py --peptides peptides.tsv --out library.tsv --n-frag 6` builds the TSV from
a tab-separated peptide list with columns `sequence, charge, protein, irt`. Check the peptide count after conversion:

```bash
TargetedFileConverter -in library.tsv -in_type tsv -out library.TraML -out_type TraML
grep -c "<Peptide " library.TraML   # must equal the number of distinct sequence+charge precursors
OpenSwathDecoyGenerator -in library.TraML -out library_decoy.TraML -method pseudo-reverse
```

The default `-method shuffle` has no seed flag and is NOT reproducible: two runs on identical
input produce different decoy peptide sequences. Use `-method pseudo-reverse` when decoys must be
reproducible (5/5 repeated runs byte-identical). `-method reverse` gives reproducible decoy
sequences and fragments, but repeated runs are not always byte-identical (3 distinct hashes in 8
runs): only the last digits of the isolation-window target m/z metadata float vary. `-method shift`
is listed by `--helphelp` but rejected every peptide as a duplicate in testing (OpenMS 3.5.0)
because it leaves the amino-acid sequence unchanged; do not rely on it.
