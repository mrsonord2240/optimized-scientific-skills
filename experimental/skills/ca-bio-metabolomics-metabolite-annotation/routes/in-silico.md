# Candidate ranking

```bash
bash scripts/run_metfrag.sh MetFragCommandLine.jar peaklist.txt candidates.csv NEUTRAL_MASS sample out
```

Replace the JAR, neutral mass and input paths. Peak list is whitespace-separated m/z and intensity without a header; LocalCSV requires `Identifier,MolecularFormula,MonoisotopicMass,InChI,InChIKey,SMILES,Name`, with comma-containing fields quoted. Source settings use 0.01 Da and 10 ppm fragment tolerances; verify suitability and installed version. The script writes parameter and result files and checks for zero candidates; existing outputs may be overwritten. A score ranks provided candidates and does not resolve all isomers. SIRIUS formula, fingerprint, structure and CANOPUS commands vary by installed version and may require a license; the source chain is not a complete authenticated file-based workflow here. Inspect current CLI help before adapting.

Done when candidate list, tolerances, alternatives and evidence limitations are reported.
