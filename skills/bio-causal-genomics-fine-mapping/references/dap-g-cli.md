# DAP-G CLI

## DAP-G Fine-Mapping

**Goal:** Deterministic posterior approximation, fast at QTL throughput (whole-transcriptome eQTL/sQTL/caQTL scans) or as an independent confirmation at a single GWAS locus.

**Approach:** DAP-G accepts the same shape of summary-statistic input as FINEMAP / susie_rss (per-SNP z-score or effect estimate, plus an LD matrix), or individual-level genotype/phenotype data in its own "sbams" text format for QTL mapping. No header row in either summary-stat input file.

```bash
# Summary statistics: snp_id z_value (two columns, no header) + a headerless,
# space-separated square LD matrix in the same SNP order
dap-g -d_z locus.zval.dat -d_ld locus.LD.dat -t 4 -o locus.dap.out
```

Full pipeline building `locus.zval.dat` / `locus.LD.dat` from a GWAS locus (reuses the same `.z`/in-sample-LD construction as `references/finemap-cli.md`): `examples/dapg_finemap.sh`.

**Checked on DAP-G (`xqwen/dap`, built from source with GSL + OpenMP, no released version tag):**
- `dap-g` **exits with code 1 on a fully successful run** — judge success by stdout (`Independent association signal clusters` block present), never the exit code, matching this Skill's general "judge by output" rule.
- Verified against the package's own bundled real example (`dap_src/sample_data/sim.1.zval.dat` + `sim.1.LD.dat`, 1001 SNPs): recovers 4 independent signal clusters at cluster_pip 0.998/0.998/0.978/0.011.
- Verified against a from-scratch synthetic GWAS locus (300 SNPs, one planted causal, in-sample LD): the planted SNP is DAP-G's top non-null candidate and heads its own signal cluster, consistent with the FINEMAP and SuSiE recovery on the same locus.
- Individual-level `-d sbams_file` mode was not exercised this session (needs a real sbams-format QTL dataset, which is out of the "single GWAS locus" scope this Skill's other worked examples use); score that code path by the tool's own documented format (`dap-g` repo `README.md` section 3.1) if it comes up.

**Output columns** (`-o` file / stdout per-SNP table): index, SNP, PIP, log10-BF-like score, cluster id (`-1` = not in a credible cluster), effect-direction sign, and a model-inclusion-frequency-style column. The `Independent association signal clusters` block at the end lists each cluster's member count, `cluster_pip`, and pairwise `average_r2` — this is DAP-G's equivalent of a SuSiE credible set.

**TORUS** (hierarchical functional-enrichment priors, paired with DAP-G for QTL scans) is a separate `xqwen/torus` repository, not bundled with `dap`; not installed or verified this session — score it by documentation if a request needs it.
