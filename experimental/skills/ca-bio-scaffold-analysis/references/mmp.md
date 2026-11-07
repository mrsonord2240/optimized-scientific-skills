## Matched Molecular Pair Analysis (MMPA) via mmpdb

**Goal:** Mine a SAR dataset for substructure transformations and their associated activity changes.

**Approach:** Fragment all compounds into core + variable side; index pairs differing by one transformation; report delta(activity) per transformation.

```bash
mmpdb fragment data.smi -o data.fragments
mmpdb index data.fragments -o data.mmpdb
mmpdb loadprops -p props.tsv data.mmpdb
mmpdb transform --smiles 'COc1ccccc1' --property pIC50 data.mmpdb
```

`props.tsv` must be a tab-separated file with an `ID` column matching the compound IDs in `data.smi` and a numeric `pIC50` column, for example:

```text
ID	pIC50
cmpd-001	6.4
cmpd-002	7.1
```

`mmpdb transform` emits TSV rows in tool order, not a ranked confidence table. Sort the output explicitly for the question at hand; useful reported columns are property-prefixed (for example, `pIC50_count`, `pIC50_avg`, `pIC50_std`, `pIC50_paired_t`, and `pIC50_p_value`). There is no mmpdb `confidence` column.

Interpret transformation effects from pair count, chemical-context diversity, dependence among pairs, uncertainty intervals, and prospective validation. Do not convert a universal pair-count/effect-size table into reliability labels.

## Context-Based MMPA

Classical MMPA: "Me -> F always +0.5 log units."
Context-based MMPA: "Me -> F adjacent to amide is +0.5; Me -> F adjacent to ester is -0.1."

Matched-pair effects can depend strongly on the local chemical environment, so report the transformation together with its attachment-point context rather than treating a global mean as universal (Raut & Dixit 2025). Use mmpdb's stored environments or a custom stratified analysis to compare context-specific effects.
