Moved out of `SKILL.md`: read when running MAGMA gene-based or gene-set analysis, or when MAGMA output feeds PoPS.

## MAGMA Gene-Based and Gene-Set Pipeline

**Goal:** Compute gene-level p-values from GWAS summary statistics and test gene sets (e.g. MSigDB pathways) for enrichment.

**Approach:** Pre-format the SNP-to-gene annotation (per-gene SNP membership using a configurable window); run gene-based analysis with `--gene-results`; downstream, test gene sets via `--set-annot`. MAGMA's lambda-correction handles LD via the reference panel.

The three steps (annotate, gene-based, gene-set) are in `scripts/magma_genebased.sh`; edit the variables at its top (GWAS file, gene-loc, reference bfile, window) and run:

```bash
bash scripts/magma_genebased.sh   # Step 1 annotate window=35,10 -> Step 2 gene-based -> Step 3 --set-annot (skipped below 200 genes)
```

The `--gene-annot` window choice is the dominant methodological lever; 35kb upstream + 10kb downstream is the FUMA recommendation but is not universally accepted. Sensitivity over 0+0, 35+10, and 50+50 is good practice for high-stakes reports.

**Step 3 has an undocumented minimum gene count.** MAGMA's `--set-annot` competitive regression conditions on 6 internal covariates (gene size, log(gene size), gene density, log(gene density), inverse MAC, log(inverse MAC)). With too few genes in `--gene-results` relative to those 6 covariates, the regression is structurally unidentifiable: MAGMA aborts with `ERROR: insufficient degrees of freedom to run analyses` (verified on a real 5-gene locus run), or, with even fewer genes, fails earlier with `no input variables to analyse` from a gene-set-variance check (verified on a real 3-gene run). Rule of thumb: gene-set enrichment needs several hundred+ genes -- **a single-locus MAGMA run (a handful of genes) should skip Step 3 entirely**; gene-set enrichment is a genome-wide- or many-loci-scale analysis, not a per-locus one. `scripts/magma_genebased.sh` checks the gene count and skips Step 3 automatically below this threshold. Windows note: this MAGMA 1.10 build also writes `<prefix>.genes.out.txt` (extra `.txt`) instead of `<prefix>.genes.out`; `scripts/magma_genebased.sh` resolves this automatically, and the same mismatch matters for PoPS below.
