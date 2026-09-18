"""
Generate the tiny offline smoke-test fixture shipped alongside this script
(toy_gwas_sumstats.tsv, toy_gene_loc.txt, toy_snp_loc.txt).

Real SNP identities and positions (chr1:55,400,000-55,600,000, hg19/GRCh37) are drawn
from the actual 1000 Genomes EUR reference .bim used by MAGMA (g1000_eur), so the
fixture is usable directly against any real 1000G EUR PLINK reference -- only the
GWAS p-values are synthetic. The window spans the real PCSK9 locus
(chr1:55,505,221-55,530,525) with one decoy gene bin upstream and one downstream,
and plants a strong synthetic signal in PCSK9 so the "true" effector gene is known
in advance. This mirrors the construction this Skill's own audit used
(GPTomics/bioSkills audit run, 2026-09-17), scaled down to ~1,800 SNPs / 3 genes so
it is a fast, tiny sanity check rather than a real analysis.

Run once with a real g1000_eur.bim (any chromosome-1-inclusive PLINK .bim; adjust
BIM_PATH) to regenerate; the committed output files in this directory do not need
regeneration to be used as-is.
"""
import numpy as np
import pandas as pd
from scipy.stats import norm

BIM_PATH = "g1000_eur.bim"  # point at your local 1000G EUR .bim (genome-wide or chr1-only)

rng = np.random.default_rng(42)

bim = pd.read_csv(
    BIM_PATH, sep=r"\s+", header=None,
    names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
)
bim = bim[(bim.CHR == 1) & (bim.BP >= 55400000) & (bim.BP <= 55600000)].reset_index(drop=True)
print(f"Loaded {len(bim)} real reference SNPs in chr1:55.4-55.6Mb window")

genes = pd.DataFrame([
    {"GENE": "DECOY_UPSTREAM",   "CHR": 1, "START": 55400000, "STOP": 55505220},
    {"GENE": "PCSK9",            "CHR": 1, "START": 55505221, "STOP": 55530525},  # real coords; planted true effector
    {"GENE": "DECOY_DOWNSTREAM", "CHR": 1, "START": 55530526, "STOP": 55600000},
])
genes["STRAND"] = "+"
genes.to_csv(
    "toy_gene_loc.txt", sep="\t", header=False, index=False,
    columns=["GENE", "CHR", "START", "STOP", "STRAND", "GENE"],
)
print("Wrote toy_gene_loc.txt (3 gene bins, MAGMA --gene-loc format)")

# Planted signal window: PCSK9 body +/- 35kb upstream / 10kb downstream (FUMA MAGMA window)
sig_lo, sig_hi = 55505221 - 35000, 55530525 + 10000
in_signal = (bim.BP >= sig_lo) & (bim.BP <= sig_hi)
print(f"SNPs inside planted PCSK9 signal window: {in_signal.sum()} / {len(bim)}")

p = np.empty(len(bim))
p[~in_signal.values] = rng.uniform(1e-4, 1.0, size=(~in_signal).sum())
z_signal = rng.normal(loc=6.5, scale=1.2, size=in_signal.sum())
p_signal = 2 * norm.sf(np.abs(z_signal))
p[in_signal.values] = np.clip(p_signal, 1e-300, 1.0)

sumstats = pd.DataFrame({
    "SNP": bim.SNP, "CHR": bim.CHR, "BP": bim.BP,
    "A1": bim.A1, "A2": bim.A2, "P": p, "N": 100000,
})
sumstats.to_csv("toy_gwas_sumstats.tsv", sep="\t", index=False)
print("Wrote toy_gwas_sumstats.tsv:", sumstats.shape)

sumstats[["SNP", "CHR", "BP"]].to_csv("toy_snp_loc.txt", sep="\t", header=False, index=False)
print("Wrote toy_snp_loc.txt")
