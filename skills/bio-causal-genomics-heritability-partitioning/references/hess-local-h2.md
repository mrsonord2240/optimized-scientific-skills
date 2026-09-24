# HESS local heritability (reference for heritability-partitioning SKILL.md)

## HESS Local h2 Pipeline

**Goal:** Estimate per-locus heritability genome-wide using LDetect partition.

**Approach:** Step 1 computes quadratic forms per locus from LD reference; step 2 estimates h2; output per-locus h2 with SE.

```bash
# HESS step 1: per-chromosome local h2 quadratic forms (run per chromosome)
for chr in {1..22}; do
    hess.py \
        --local-hsqg trait.sumstats.gz \
        --chrom $chr \
        --bfile 1000G_EUR_chr${chr} \
        --partition LDetect_EUR_chr${chr}.bed \
        --out hess_chr${chr}
done

# HESS step 2: aggregate across chromosomes and estimate h2
hess.py \
    --prefix hess_chr \
    --out trait_local_h2 \
    --tot-hsqg <total_h2_estimate> <total_h2_SE>
# Provide total h2 and SE from LDSC for the global constraint
```

LDetect partition files (Berisa & Pickrell 2016) are available pre-computed for EUR/EAS/AFR at https://bitbucket.org/nygcresearch/ldetect-data. HESS detects high-h2 loci suitable for fine-mapping prioritization.

Sumstats file needs columns `SNP CHR BP A1 A2 Z N` (tab- or whitespace-delimited, single-letter
alleles, rsIDs); `--chrom` must match the `CHR` value in the sumstats file exactly (string match, not
numeric). Run step 1 for one chromosome as
`bash examples/hess_local_h2.sh <sumstats.txt> <chrom> <bfile_prefix> <partition.bed> <out_prefix>`;
loop it per chromosome, then run step 2 (`hess.py --prefix ... --tot-hsqg <h2> <SE> --out ...`)
yourself once every chromosome's step 1 has finished -- step 2 requires all 22 autosomes' step-1
output files to be present, by design.

## Failure Mode: HESS locus instability at sparse SNP density

**Trigger:** Running HESS at a locus with < 1000 LD-pruned SNPs (e.g. centromere-adjacent region).

**Mechanism:** HESS quadratic form on projected effect estimates is unstable at low SNP density; matrix conditioning explodes.

**Symptom:** Locus h2 estimate negative or > 0.5 (unphysical); standard error very large; subsequent loci stable.

**Fix:** Require >= 1000 SNPs per locus; use LDetect partition (Berisa & Pickrell 2016 Bioinformatics 32:283) which targets ~1700 loci genome-wide; discard sparse loci or merge with neighbors.

## Install

```bash
# HESS 0.5.4-beta (huwenboshi/hess). Checked 2026-09-21: there is no separate "Python 3 branch" --
# the repo has one branch (hess-0.5, HEAD c7a2330) and no requirements.txt. The shipped code does not
# run under Python 3 as-is; five mechanical patches make it run (verified end-to-end on real genotypes
# and real per-SNP regression Z-scores; see examples/hess_local_h2.sh):
git clone https://github.com/huwenboshi/hess.git
pip install numpy pandas scipy pysnptools   # pandas 3.x is fine once patch 4 below is applied
cd hess
# 1. src/estimation.py: implicit relative imports fail under Python 3
sed -i "s/^from sumstats import \*/from .sumstats import */" src/estimation.py
sed -i "s/^from refpanel import \*/from .refpanel import */" src/estimation.py
# 2. np.float was removed in numpy>=1.24 (src/estimation.py, 4 occurrences)
sed -i 's/np\.float)/float)/g' src/estimation.py
# 3. xrange is Python-2-only (src/estimation.py, src/refpanel.py, src/sumstats.py)
sed -i 's/\bxrange\b/range/g' src/estimation.py src/refpanel.py src/sumstats.py
# 4. pandas>=2.2 removed delim_whitespace (src/estimation.py, src/refpanel.py)
#    replace each `delim_whitespace=True` with `sep=r"\s+"`
# 5. gzip.open(..., 'w'/'r') defaults to binary mode in Python 3 but the code writes/reads str
#    (src/estimation.py): change every gzip.open(..., 'w') to 'wt' and every (..., 'r') to 'rt'
```

Verified 2026-09-21 with these five patches applied, on a synthetic locus (real 957-individual,
1129-SNP chr1 panel, 20 planted causal SNPs, real per-SNP OLS Z-scores): step 1
(`--local-hsqg`) correctly reports 505 and 514 SNPs in the two halves of the locus and writes real
`.info.gz`/`.eig.gz`/`.prjsq.gz` output -- see `examples/hess_local_h2.sh`. Step 2 was not exercised
end-to-end in that session (it requires step-1 output for all 22 autosomes, by the tool's own design,
which a single-chromosome synthetic test cannot supply); read `local_hsqg_step2_helper` in
`src/estimation.py` before relying on step 2's output on real genome-wide data.
