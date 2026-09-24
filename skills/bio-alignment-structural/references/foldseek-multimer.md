# Foldseek-Multimer: complex search and clustering

### Foldseek-Multimer for Database-Scale Complex Search

For multi-chain complex search over thousands of complexes or more (a folder of complexes, or a `foldseek createdb` of PDB entries), US-align is too slow; Foldseek-Multimer (Kim et al 2025 Nat Methods) is the modern default. `foldseek databases` lists no multimer database, so the target is a folder or a database you build. Two modes share the chain-pairing prefilter:

| Mode | Algorithm | Use when |
|------|-----------|----------|
| `Foldseek-MM` (default) | 3Di+AA Gotoh per chain pair | Fast database search; default speed-sensitivity tradeoff |
| `Foldseek-MM-TM` | TM-align per chain pair after prefilter (`--alignment-type 1`) | Top-hit refinement with full TM-score |

```bash
foldseek easy-multimersearch query_complex.pdb complexes_dir/ result tmp/
foldseek easy-multimersearch query_complex.pdb complexes_dir/ result tmp/ --alignment-type 1   # Foldseek-MM-TM

# Cluster: pass ONE DIRECTORY of complexes. A file list or *.pdb clusters only the last file and still exits 0.
foldseek easy-multimercluster complexes_dir/ cluster_result tmp/ --multimer-tm-threshold 0.65
```

`result` holds chain-level hits. The complex-level scores are in `result_report` (tab-separated: query complex, target complex, query chains, target chains, qTM, tTM, rotation `u`, translation `t`, assembly id); rank and filter on columns 5-6, since `easy-multimersearch` rejects `--multimer-tm-threshold` (that flag exists only on `easy-multimercluster`) and its `--tmscore-threshold` filters chain alignments only. As under `--alignment-type 1` elsewhere, E-values in `result` are not meaningful there. Check: searching 1IRD (dimer) against a directory holding the 1A3N tetramer gives qTM 0.981 / tTM 0.492 for chains A,B vs A,B, matching US-align (0.977 / 0.494).

Decision guide: pairwise complex pair on a few hundred targets -> US-align with `-mm 1 -ter 0`. Database search across thousands to millions of complex entries -> Foldseek-Multimer: the paper reports 3-4 orders of magnitude over US-align at that scale. It is not faster by that factor on small sets: in one audit run on 62 small oligomers it was 1.7x (folder, including `createdb`) to 5.8x (prebuilt database) faster than US-align, and its prefilter reported only 22 of the 62 targets, though no pair with TM > 0.5 was missing.

**Reporting convention:** Always quote BOTH the multimer TM-score (whole-complex) AND the worst per-chain TM-score; high multimer TM with one low per-chain score signals topology-matched but locally divergent chains (often promiscuous binders or paralog swaps), which is biologically distinct from a uniformly-conserved complex.
