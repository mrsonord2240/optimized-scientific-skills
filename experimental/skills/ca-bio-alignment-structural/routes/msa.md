# Structural multiple alignment

```bash
foldmason easy-msa structures/*.pdb result tmp/ --refine-iters 100 --refine-seed 42 --report-mode 1
```

Use at least two structures. Rows represent chains, not files. Set a refinement seed or use zero iterations for reproducible results. The amino acid and 3Di MSAs, guide tree and confidence report are described in `references/foldmason.md`; MUSTANG is an alternative in `references/search.md`. Downstream MSA statistics use `bio-alignment-msa-statistics`.

Done when the requested MSA and confidence output are produced and chain counts are checked.
