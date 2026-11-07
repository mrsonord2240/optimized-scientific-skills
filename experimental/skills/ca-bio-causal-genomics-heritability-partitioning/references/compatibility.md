# Tool Install Notes
```bash
# LDSC: use CBIIT/ldsc (NIH CBIIT's maintained Python 3 fork; bulik/ldsc's own README
# points here as of 2026-01-16). Checked on commit 1f09cf0 (2026-08-10), Python 3.9.
git clone https://github.com/CBIIT/ldsc.git
cd ldsc
conda create --name ldsc python=3.9 -y && conda activate ldsc   # or: micromamba create -n ldsc -c conda-forge -c bioconda python=3.9 bitarray=2 pybedtools=0.10.0 -y
pip install numpy==1.21.5 pandas==1.3.3 scipy==1.7.3   # environment3.yml's pins; conda-installed numpy/pandas float otherwise
./ldsc.py -h   # verify: prints the full flag list, no traceback

# --h2 and --rg run as-is against this clone. --h2-cts needs one line patched first:
#   ldscore/sumstats.py ~line 285: ref_ld_cts = ... .loc[:,1:])  ->  ... .iloc[:,1:])
# (a leftover from replacing pandas' removed .ix[:,1:] with the label-based .loc instead of
# the positional .iloc; pandas >=1.0 raises TypeError on the unpatched line.) Without the
# patch, --h2-cts fails with "cannot do slice indexing on Index with these indexers [1] of
# type int". LD score files must be gzipped (`.l2.ldscore.gz`) -- ldscore() hardcodes that
# suffix; every official reference bundle already ships gzipped.

# Pre-computed reference resources (one-time download)
# EUR LD scores (HapMap3 SNPs)
wget https://alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr.tar.bz2
# Baseline-LD v2.2 (Gazal 2017)
wget https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
# 1000G EUR frequency files
wget https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_frq.tgz
# Weights
wget https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz
# Multi-tissue chromatin ldcts (Finucane 2018)
wget https://alkesgroup.broadinstitute.org/LDSCORE/Multi_tissue_chromatin_1000Gv3_ldscores.tgz
```

Before pointing LDSC at real GWAS data, confirm the install works: `LDSC_DIR=./ldsc bash the source checkout smoke-test script`
runs `--h2`, `--rg`, and `--h2-cts` against the simulated fixtures in the clone's own `test/` directory
and prints real, checkable regression output for each.

LDAK, HESS and HDL install notes live in their `references/` files (see "Per-Method Reference Files").

```bash
# BOLT-LMM / BOLT-REML (checked 2026-09-21: the old v2.4.1 URL below 404s; current release is v2.5,
# older versions moved under .../downloads/old/)
wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.5.tar.gz
# Segfaults immediately (even --help) on this Skill's own verification environment (WSL, 2026-09-21) --
# not yet root-caused; verify with `bolt --help` before trusting the binary on a new machine.
```

GCTA-GREML and Popcorn each have their own pipeline, install and worked example: see "Per-Method
Reference Files" below (`references/gcta-greml.md`, `references/popcorn-transancestry.md`).
