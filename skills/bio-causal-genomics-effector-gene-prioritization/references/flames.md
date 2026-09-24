Moved out of `SKILL.md`: read when running FLAMES on a fine-mapped credible set.

## FLAMES

**Goal:** Score candidate genes at a fine-mapped locus with FLAMES' pretrained XGBoost classifier (multimodal V2G evidence: VEP, CADD, eQTL colocalization, Hi-C, enhancer-gene links) combined with PoPS, producing a calibrated precision estimate per gene.

**Install** (checked 2026-09-21, `Marijn-Schipper/FLAMES` HEAD `159e83a`, version 1.1.3):

```bash
git clone https://github.com/Marijn-Schipper/FLAMES
cd FLAMES
conda env create -f environment.yml   # python 3.8.13, scikit-learn 1.2.2, xgboost 1.7.5, pandas 2.0.3
conda activate FLAMES
pip install -r requirements.txt

# Annotation data (~1.7 GB) -- public, unauthenticated
curl -L -o Annotation_data.tar.gz "https://zenodo.org/records/12635505/files/Annotation_data.tar.gz?download=1"
tar xzf Annotation_data.tar.gz   # -> Annotation_data/{ENSG,eQTLGen,GTEx_v8_finemapping_CAVIAR,eQTL_catalog,
                                  #                     HACER,Cicero_whole_blood,Jung_PCHiC,Javierre_PCHiC,
                                  #                     ABC_EP,ABC_CRISPR,EP_correlation_F5,GeneHancer,
                                  #                     RoadMapEpi,EpiMap,Promoter_regions}
```

**Run** (real end-to-end run, verified 2026-09-21, on FLAMES' own bundled `example_data/` -- the four dizygotic-twinning GWAS loci from the FLAMES paper; MAGMA gene Z, MAGMA tissue enrichment, PoPS preds, and per-locus 95%-pruned SuSiE credible sets all ship with the repo):

```bash
cd FLAMES/example_data
python ../FLAMES.py annotate \
    -a ../Annotation_data/ -p PoPS.preds -m magma.genes.out \
    -mt magma_exp_gtex_v8_ts_avg_log2TPM.txt.gsa.out -id indexfile.txt \
    -pc prob1 -sc cred1 -g genes.txt -c95 False
# writes annots/annotated_locus_{1..4}.txt (about 20-80s/locus; queries the public VEP + CADD REST APIs
# per credible-set variant unless -cv/-vc point at a local VEP install)

python ../FLAMES.py FLAMES -id indexfile.txt -o ./
# writes FLAMES_scores.pred (genes above the calibrated 75% cumulative-precision threshold) and
# FLAMES_scores.raw (every scored gene)
```

Real result (`FLAMES_scores.pred`): the correct twinning-associated gene at each of the 4 loci --
`GNRH1` (locus 1, precision 0.88), `FSHB` (locus 2, precision 0.87), `SMAD3` (locus 3, precision 0.95),
`ZFPM1` (locus 4, precision 0.98) -- each the reproductive-hormone or TGF-beta-pathway gene the FLAMES
paper itself nominates at these loci, not merely the nearest gene (`FLAMES_scores.raw` shows several
non-causal candidate genes per locus scoring lower).

**On your own data**, you need: MAGMA `.genes.out`/`.genes.raw` (`references/magma-gene-based.md`),
MAGMA tissue-expression `.gsa.out` (`--gene-covar` against the bundled GTEx file from the same Zenodo
record), PoPS `.preds` (`references/pops.md`), and credible sets as two tab-separated columns
(`CHR:BP:A1:A2` or `CHR:BP:A1_A2`, then PIP) -- format matches this Skill's other credible-set inputs
(cross-reference causal-genomics/fine-mapping). Genome build defaults to hg19/GRCh37; pass `-b GRCh38`
for hg38 credible sets (FLAMES lifts over internally via `pyliftover`).

**Speed / API note:** default `annotate` queries the public Ensembl VEP and CADD REST APIs once per
credible-set variant -- this is the bottleneck (verified: ~10-80s/locus on a 1-20 variant credible set,
scaling with credible-set size and network latency). For larger batches, pass `--cmd-vep`/`--vep-cache`
(command-line VEP) and `-t`/`-cf` (tabixed local CADD) to avoid the live API calls entirely.
