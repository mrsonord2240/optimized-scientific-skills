# Foldmason easy-msa and per-column LDDT

### Foldmason easy-msa

```bash
foldmason easy-msa structures/*.pdb result tmp/ \
    --refine-iters 100 --refine-seed 42 \
    --report-mode 1

# Outputs: result_aa.fa (amino-acid MSA), result_3di.fa (3Di MSA), result.nw (guide tree), result.html (LDDT report)
```

`--refine-iters` controls iterative MSA refinement (default 0, which is deterministic). Refinement is random unless `--refine-seed` is set: unseeded `--refine-iters 100` runs give different MSAs (how different depends on the set: from under 1% to about 40% of aligned pairs differing between runs), while `--refine-seed 42` gave byte-identical FASTA on repeat. Set a seed or use `--refine-iters 0` for any result you compare or report. `--report-mode 1` produces an HTML report with per-column LDDT confidence. The 3Di MSA can be used directly for evolutionary analyses where structure rather than sequence is the appropriate signal. Foldmason writes one row per chain, named `<file>_<chain>` (a single-chain file keeps its bare stem): 11 files gave 18 rows, so count rows as chains, not structures.

**Per-column LDDT extraction:** Run with `--report-mode 2` to produce machine-readable JSON output:

```bash
foldmason easy-msa structures/*.pdb result tmp/ --report-mode 2
# Produces result.json: keys entries (rows: name, aa, ss, ca), scores, tree, statistics
```

```python
import json
with open('result.json') as f:
    report = json.load(f)
lddt_per_column = report['scores']          # NOT 'per_column_lddt' (that key does not exist; .get returns None)
assert len(lddt_per_column) == len(report['entries'][0]['aa'])   # one value per MSA column
```

`-1` marks a column with no LDDT (113 of 486 columns in an 11-structure globin/kinase run). Checked on Foldmason 4.dd3c235; inspect a sample `result.json` if the version differs.
