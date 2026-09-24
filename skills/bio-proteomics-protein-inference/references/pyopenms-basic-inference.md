# Group Proteins with pyOpenMS: Basic Inference + Greedy Resolution + Picked-Group FDR

### Group Proteins with pyOpenMS (aggregation + greedy resolution)

**Goal:** Turn an FDR-filtered peptide identification list into protein groups with a leading protein, resolving shared-peptide ambiguity.

**Approach:** Load the idXML from peptide identification, run `BasicProteinInferenceAlgorithm` (score aggregation per protein) with indistinguishable-group annotation AND greedy group resolution on -- without resolution, subsumable and shared-only proteins stay as their own groups -- apply picked protein-group FDR with the built-in, and read the groups off the protein identification run as records: `leading_protein`, `accessions`, `n_peptides`, `n_unique_peptides` (unique to the GROUP), `is_decoy`, `qvalue` -- the contract to bind downstream code to, and the same record `examples/protein_groups.py` returns.

```bash
python scripts/group_proteins_pyopenms.py peptides_1pct_fdr.idXML --decoy-prefix DECOY_ --out groups.tsv
```

`scripts/group_proteins_pyopenms.py` sets `annotate_indistinguishable_groups` and `greedy_group_resolution`, applies `applyPickedProteinFDR`, builds the group records, and prints the groups passing `--fdr` (default 0.01). `--decoy-prefix` defaults to `DECOY_`; set it to the search's prefix: `DECOY_` OpenMS/Comet, `rev_` Philosopher, `REV__` MaxQuant.
