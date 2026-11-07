# Per-Tool Failure Modes

### LD-induced TWAS false positives (most common pitfall)

**Trigger:** Two or more genes at the same locus have correlated cis-eQTLs (shared causal eQTL SNP or LD-linked eQTL SNPs).

**Mechanism:** TWAS Z-scores are linear combinations of SNP Z-scores weighted by per-gene SNP effects. When two genes share many high-weight SNPs (e.g. nearby genes regulated by the same enhancer or LD-tagged independent eQTLs), their TWAS Z-scores are positively correlated. A single causal GWAS variant therefore produces significant Z at multiple co-regulated genes (Wainberg 2019 Nat Genet 51:592; Mancuso 2019 Nat Genet 51:675).

**Symptom:** A GWAS lead locus shows 3-10 genes all passing genome-wide TWAS significance (p < 2.3e-6 ~ 0.05/22k); per-gene LocusZoom-style plots look near-identical; conditional analysis (`FUSION.post_process.R`, see `fusion-pipeline.md`) reveals only 1-2 independent gene signals; the genes lie within 1 Mb of each other.

**Fix:** Always run FOCUS after TWAS to obtain per-gene PIPs. Report only genes with PIP >= 0.8 as candidate causal; report co-significant genes with PIP < 0.5 as LD-tagged. Cross-check with cis-eQTL coloc (PP.H4 >= 0.7) for the candidate causal gene.

### Tissue mis-specification

**Trigger:** Running TWAS in a tissue that does not host the causal regulatory effect (e.g. whole blood for a psychiatric trait; pancreas for an LDL trait).

**Mechanism:** A gene's cis-eQTL effect size varies across tissues; a wrong-tissue model has weaker per-gene prediction R^2 and lower power. Conversely, eQTL effects in the wrong tissue can still tag the GWAS signal via LD and produce spurious associations not present in the causal tissue.

**Symptom:** Strong TWAS signal in a tissue biologically irrelevant to the trait; null in the expected tissue; tissue-prioritisation methods (LDSC-SEG, Finucane 2018 Nat Genet 50:621; RolyPoly (Calderon 2017 AJHG 101:686)) disagree with the TWAS tissue.

**Fix:** Run S-MultiXcan to combine tissues if causal tissue is unknown. For prioritisation, use LDSC-SEG / CELLEX / EWCE on the GWAS sumstats independently of TWAS, and report TWAS in the prioritised tissues. Never report a single-tissue TWAS hit as causal without independent tissue evidence (single-cell eQTL, chromatin accessibility in matched cell type).

### Ancestry mismatch in prediction weights

**Trigger:** Running TWAS on a non-EUR GWAS using GTEx (~ 85% EUR) weights, or vice versa.

**Mechanism:** Cis-eQTL effect sizes and LD structure are ancestry-specific; prediction weights trained in one ancestry transfer with reduced R^2 and biased Z-scores in another. Power is lost preferentially at loci where the causal eQTL is not shared across ancestries (Patel 2022 AJHG 109:1286).

**Symptom:** Genome-wide TWAS hit count much lower than expected given GWAS power; non-EUR-specific GWAS loci fail to produce TWAS hits; per-gene prediction R^2 substantially reduced.

**Fix:** Use ancestry-matched prediction panels where available: MESA multi-ethnic eQTL (Mogil 2018 PLoS Genet), eQTLGen-Asian, AFGR (Africa) when published, or MAGE (Taliun-style multi-ancestry eQTL). Move to MA-FOCUS for cross-ancestry joint fine-mapping. Document the ancestry assumption explicitly in methods.

| Panel | N donors | Tissues / cells | Ancestry | Use case |
|-------|----------|-----------------|----------|----------|
| GTEx v8 (full) | 838 | 49 tissues | EUR (~85%) | Default standard for general TWAS |
| GTEx v8 MASHR-EUR | -- | 49 tissues | EUR | Primary; sparser SNP set per gene |
| eQTLGen | 31,684 | Whole blood | EUR (>95%) | Highest blood power; cis + trans available |
| MESA Monocytes (Mogil 2018) | 1,163 | Monocytes | Multi-ancestry | AFR / HIS-relevant analyses |
| MESA Monocytes-AFA | 233 | Monocytes | AFR | AFR-specific immune traits |
| AFGR | ~2,000 | Whole blood | AFR | AFR (emerging; release-dependent) |
| OneK1K (Yazar 2022) | 982 | PBMC, 14 cell types | EUR | sc-TWAS in immune cell types |
| PsychENCODE (Wang 2018) | ~1,300 | Prefrontal cortex | EUR | Neuropsychiatric traits |
| BrainSeq Phase 2 | ~350 | DLPFC, hippocampus | EUR + AFR | Neuropsychiatric replication |

### Low-N tissue weights are unstable

**Trigger:** Using a GTEx tissue with N < 100 donors (e.g. several brain sub-regions, kidney cortex in v7).

**Mechanism:** Per-gene elastic-net weights are cross-validated with the available donors. Below ~ 100 donors, the cross-validation R^2 has high variance and the heritability filter (FUSION requires hsq_p < 0.01) drops many genes. Surviving weights overfit, inflating per-gene Z under the null.

**Symptom:** Tissue produces unusually high TWAS hit count or unusually high genomic inflation; per-gene CV R^2 distribution is bimodal with a long heavy tail.

**Fix:** Skip GTEx tissues with N < 100 unless biologically essential. Substitute eQTLGen for whole blood (N ~ 31k, Vosa 2021 Nat Genet 53:1300) where blood is acceptable. For brain, use PsychENCODE (N ~ 1300, Wang 2018 Science 362:eaat8464) or BrainSeq (N ~ 350+) when available; verify the matching prediction-weight panel exists.

GTEx v8 small-N tissues (skip those below 100 unless biologically required; the two brain rows are marginal):

| Tissue | GTEx v8 N |
|--------|-----------|
| Kidney Medulla | 4 |
| Cervix - Endocervix | 10 |
| Cervix - Ectocervix | 9 |
| Fallopian Tube | 9 |
| Bladder | 21 |
| Brain - Substantia nigra | 139 |
| Brain - Spinal cord (cervical c-1) | 159 |

### HLA region

**Trigger:** Any gene within chr6:25-35 Mb (hg38; extended MHC) reported by TWAS.

**Mechanism:** Long-range LD (r2 > 0.5 over many Mb) and extreme structural variation mean per-gene prediction weights at HLA capture haplotype rather than gene-specific regulation. Standard TWAS gene-level inference is biologically meaningless here.

**Fix:** Exclude chr6:25-35 Mb from genome-wide TWAS summaries by default. For HLA-driven traits (autoimmune, infection, transplantation), impute classical HLA alleles using one of:

| Tool | Reference | Notes |
|------|-----------|-------|
| HIBAG | Zheng 2014 Pharmacogenomics J 14:192 | R package; pre-trained per-ancestry classifiers |
| SNP2HLA | Jia 2013 PLoS One 8:e64683 | Beagle-based imputation; supports T1DGC reference |
| HLA-TAPAS | Luo 2021 Nat Genet 53:1504 | Current standard; multi-ancestry reference; recommended for new analyses |

Then test classical alleles plus amino-acid residues (Raychaudhuri 2012 Nat Genet 44:291 set the gold standard for residue-level association in MHC). Do NOT run SNP-level TWAS inside the MHC.

### Correlated-expression confounding (co-regulated genes)

**Trigger:** Two or more genes are functionally co-regulated by a single TF or enhancer, producing nearly-identical predicted-expression vectors.

**Mechanism:** Even with separate per-gene cis-eQTL prediction, downstream co-regulation makes predicted expression highly correlated; TWAS cannot distinguish which gene mediates the trait.

**Symptom:** FOCUS credible gene set contains multiple genes with PIP roughly equal (e.g. three genes at 0.3 each); functional follow-up (MPRA, CRISPRi screens) is needed to break the tie.

**Fix:** Acknowledge the limit of statistical resolution; report the full credible gene set and prioritise on orthogonal evidence (CRISPRi/CRISPRa effect size in matched cell type, e.g. Open Targets-style, or MPRA at allelic series; protein-level pQTL coloc if available).
