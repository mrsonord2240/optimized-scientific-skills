# Picked Protein-Group FDR

### Picked Protein-Group FDR

**Goal:** Estimate protein-group FDR without the inflation that the reused PSM formula causes on large data.

**Approach:** For each target group, find its decoy counterpart (same accessions with the decoy prefix); keep only the higher-scoring member of each target/decoy PAIR; rank the picked set and count decoys as the FDR estimate. For idXML, prefer the built-in `FalseDiscoveryRate().applyPickedProteinFDR(prot_id, String(prefix), True, True)` shown in `pyopenms-basic-inference.md`; for group-level work on large data, the kusterlab `picked_group_fdr` package implements The 2022. The script pairs by the exact accession set, so decoy groups whose membership differs from their target's stay unpaired and are counted unpicked. The decoy prefix is tool-specific (`DECOY_` OpenMS/Comet, `rev_` FragPipe/Philosopher's `--tag` default, `REV__` MaxQuant) and must be passed explicitly (`DECOY_PREFIX` in the blocks). A wrong prefix does not pass silently in either implementation: `applyPickedProteinFDR` raises `IndexError: invalid unordered_map<K, T> key`, and the script raises `ValueError: no decoy groups with prefix ...`. Treat both as "the prefix is wrong", not as a corrupt input file.

```bash
python scripts/picked_group_fdr.py groups.tsv --decoy-prefix DECOY_ --score-col probability --out passing.tsv
```

`scripts/picked_group_fdr.py` is also importable (`from picked_group_fdr import picked_group_fdr`); it raises on a wrong prefix and warns below 10 decoy groups. The self-contained demo of parsimony plus the same pairing is `examples/protein_groups.py`.
