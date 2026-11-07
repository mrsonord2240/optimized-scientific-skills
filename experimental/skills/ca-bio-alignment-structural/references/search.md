# Search



Run Foldseek for AlphaFoldDB-scale structural search; it indexes a 20-letter 3Di alphabet for thousand- to million-fold speedup over TM-align at comparable same-fold sensitivity.

```bash
# Search query structures against AFDB (default: --alignment-type 2)
foldseek easy-search query.pdb afdb_database result.m8 tmp/

# Refine top hits with full TM-align rotation (slower but global TM-score); filter on alntmscore, not E-value
foldseek easy-search query.pdb afdb_database result.m8 tmp/ --alignment-type 1

# All-versus-all clustering at TM > 0.5
foldseek easy-cluster structures/*.pdb cluster_result tmp/ --tmscore-threshold 0.5

# Custom output columns
foldseek easy-search query.pdb afdb_database result.m8 tmp/ \
    --format-output query,target,evalue,alntmscore,qtmscore,ttmscore,lddt,bits
```

| Foldseek `--alignment-type` | Algorithm | Use when |
|-----------------------------|-----------|----------|
| 0 | 3Di Gotoh-Smith-Waterman (local) | Not recommended; 3Di alone, no amino-acid signal |
| 1 | TMalign (global) | Refine top hits with full TM-score; slow. E-values are meaningless (0.88-0.99 for true homologs, alnTM 0.93): an `evalue < 1e-3` filter returns 0 hits, so filter on `alntmscore` (and `lddt`) |
| 2 | 3Di+AA Gotoh-Smith-Waterman (local) | Default; best speed-sensitivity tradeoff |

`--max-seqs` (default 1000; `examples/foldseek_search.py` uses 200) silently caps the rows per query: a printed count equal to the cap is a cap, not the hit total.

### pLDDT-Filtering AlphaFold Structures Before Foldseek

Mask residues with pLDDT < 70 before Foldseek indexing or search; their backbone coordinates encode as random 3Di letters and contaminate hits below TM ~ 0.4. AFDB clusters from Barrio-Hernandez et al 2023 are pre-filtered; ESMFold predictions need the same step.

```bash
# Mask low-pLDDT residues during database creation (B-factor column stores pLDDT)
# The *.pdb glob requires shell expansion; for Python subprocess calls use glob.glob()
# to expand the file list before passing as argv.
foldseek createdb --mask-bfactor-threshold 70.0 *.pdb afdb_masked
```

## Structural Multiple Sequence Alignment

| Tool | Reference | Best for |
|------|-----------|----------|
| Foldmason | Gilchrist et al 2026 Science 391:485 | Billion-protein-scale structural MSA on AFDB; the structural counterpart to MAFFT |
| MUSTANG | Konagurthu et al 2006 Proteins | Pure structural MSA via residue-residue equivalences; `mustang-3.2.3 -i a.pdb b.pdb c.pdb -o out -F fasta` writes `out.afasta` + `out.pdb` (3 myoglobins: 3-row FASTA) |

Not covered, because no runnable path could be verified: T-Coffee Expresso and 3D-Coffee (T-Coffee 12.00.7: Expresso needs a BLAST server; 3D-Coffee retries retired RCSB/wwPDB endpoints and falls back to `proba_pair`), PROMALS3D and mTM-align (web servers), FATCAT, ChimeraX (licence-gated download) and the protein-language-model aligners TM-Vec, vcMSA, DEDAL and pLM-BLAST (GPU and model weights). Use Foldmason or MUSTANG.

### Foldmason easy-msa

`foldmason easy-msa structures/*.pdb result tmp/ --refine-iters 100 --refine-seed 42 --report-mode 1` writes `result_aa.fa`, `result_3di.fa`, `result.nw` and `result.html`. Always set `--refine-seed` (or `--refine-iters 0`): unseeded refinement is random. Rows are chains, not files. Per-column LDDT is under the JSON key `scores` (`--report-mode 2`), not `per_column_lddt`. Details and the extraction snippet: `references/foldmason.md`.
