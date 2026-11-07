"""Synthetic replica of variantbench `cis_mr_bcaa`.

Everything here is FABRICATED. Summary statistics are simulated from a known causal
model; positions are fictional; gene names of the BCAA catabolism pathway are real
because the task prompt names BCKDK and the pathway knowledge is part of the test.
No number here is a real-world finding.

World (design source: the original's public `notes`, not our Skills):
  * BCKDK (and BCKDHA, BCKDHB, DLD, BCAT1) have NO usable common-variant
    instruments: no cis signal in the BCAA GWAS.
  * DBT, PPM1K, BCAT2 carry strong cis BCAA signals. Their cis-MR estimates for T2D
    are small and of mixed sign, so the evidence does NOT support BCKDK inhibition.
  * Genome-wide BCAA instruments are horizontally pleiotropic: a naive genome-wide
    IVW gives a clearly positive effect (the "YES" trap).
  * T2D file uses different column names, row order and allele orientation
    (15% of rows flipped) so harmonisation is required.

Run: python generate.py  -> data/ + eval.json + truth.json
Seed 1 was chosen from seeds 1-6 because its realised cis estimates sit close to the
designed truth with no locus individually significant (a design constraint on the
fabricated data, not a search for a favourable result).
"""
import gzip
import json
import os
from math import erfc, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
SEED = int(os.environ.get("SEED", 1))
rng = np.random.default_rng(SEED)

L = 250                      # SNPs per region
N_POP, N_REF = 20_000, 500   # diploid individuals: LD population / reference panel
N_BCAA = 115_000
N_CASE, N_CTRL = 180_000, 1_100_000
PHI = N_CASE / (N_CASE + N_CTRL)
NEFF_Y = (N_CASE + N_CTRL) * PHI * (1 - PHI)
THETA_GW = 0.30              # pleiotropy-inflated genome-wide "effect"

# (name, chr, role, theta_true) -- theta = true effect of +1 SD BCAA on T2D log-odds
PATHWAY = [
    ("DBT",    1, "signal", -0.08),
    ("PPM1K",  4, "signal", +0.09),
    ("BCAT2", 19, "signal", +0.01),
    ("BCKDK", 16, "null",   0.0),
    ("BCKDHA", 19, "null",  0.0),
    ("BCKDHB", 6, "null",   0.0),
    ("DLD",    7, "null",   0.0),
    ("BCAT1", 12, "null",   0.0),
]
N_GW = 50

ALLELE_PAIRS = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"),
                ("A", "C"), ("C", "A"), ("G", "T"), ("T", "G")]   # none palindromic


def qnorm(p):
    from scipy.stats import norm
    return norm.ppf(p)


def pval(z):
    return erfc(abs(z) / sqrt(2))


def simulate_region(chrom, start):
    """Return dict with positions, alleles, freq, population R, reference genotypes."""
    pos = (start + np.cumsum(rng.uniform(800, 2400, L))).astype(int)
    maf = 0.03 + 0.47 * rng.beta(1, 1.5, L)
    p = np.where(rng.random(L) < 0.5, maf, 1 - maf)               # freq of A1
    rho = rng.uniform(0.94, 0.995, L)
    rho[rng.random(L) < 0.08] = 0.3                                # recombination hotspots
    nh = 2 * (N_POP + N_REF)
    z = rng.standard_normal(nh)
    H = np.empty((nh, L), dtype=np.int8)
    for i in range(L):
        if i:
            z = rho[i] * z + np.sqrt(1 - rho[i] ** 2) * rng.standard_normal(nh)
        H[:, i] = z < qnorm(p[i])
    G = H[0::2] + H[1::2]                                          # count of A1
    pop, ref = G[:N_POP], G[N_POP:]
    eaf = pop.mean(0) / 2
    R = np.corrcoef(pop.T.astype(np.float32)).astype(np.float64)
    pairs = [ALLELE_PAIRS[k] for k in rng.integers(0, len(ALLELE_PAIRS), L)]
    return dict(chr=chrom, pos=pos, eaf=eaf, R=R, ref=ref,
                a1=[a for a, _ in pairs], a2=[b for _, b in pairs])


def mvn_noise(R, se):
    w, V = np.linalg.eigh(R)
    w = np.clip(w, 0, None)
    return se * (V @ (np.sqrt(w) * rng.standard_normal(len(w))))


def marginal(R, sd, beta_causal):
    """per-allele marginal effects implied by sparse causal per-allele effects."""
    return (R * np.outer(1 / sd, sd)) @ beta_causal


def pick_causal(reg, centre, n, span=60_000):
    sd_ok = (reg["eaf"] > 0.15) & (reg["eaf"] < 0.85) & (np.abs(reg["pos"] - centre) < span)
    idx = np.flatnonzero(sd_ok)
    return rng.choice(idx, size=n, replace=False)


regions, genes = [], []
occupancy = {}


def next_start(chrom):
    k = occupancy.get(chrom, 0)
    occupancy[chrom] = k + 1
    return 10_000_000 + k * 8_000_000


def build(name, chrom, role, theta, is_gw=False):
    start = next_start(chrom)
    reg = simulate_region(chrom, start)
    centre = int(reg["pos"][L // 2])
    sd = np.sqrt(2 * reg["eaf"] * (1 - reg["eaf"]))
    se_x = 1 / (sd * np.sqrt(N_BCAA))
    se_y = 1 / (sd * np.sqrt(NEFF_Y))
    bx_c = np.zeros(L)
    by_c = np.zeros(L)
    if role == "signal":
        i1, i2 = pick_causal(reg, centre, 2)
        bx_c[i1] = rng.choice([-1, 1]) * rng.uniform(0.075, 0.10)
        bx_c[i2] = rng.choice([-1, 1]) * rng.uniform(0.025, 0.04)
        by_c = theta * bx_c                                        # clean cis effect, no pleiotropy
    elif role == "gw":
        for i in pick_causal(reg, centre, rng.integers(1, 3), span=150_000):
            bx_c[i] = rng.choice([-1, 1]) * rng.uniform(0.02, 0.05)
        by_c = THETA_GW * bx_c + rng.normal(0, 0.006, L) * (bx_c != 0)
    bx_m, by_m = marginal(reg["R"], sd, bx_c), marginal(reg["R"], sd, by_c)
    for attempt in range(500):
        bx = bx_m + mvn_noise(reg["R"], se_x)
        if role != "null" or np.abs(bx / se_x).max() < 2.9:        # null regions: guarantee no instrument
            break
    else:
        raise RuntimeError("could not draw a null region")
    by = by_m + mvn_noise(reg["R"], se_y)
    reg.update(name=name, role=role, theta=theta, bx=bx, se_x=se_x, by=by, se_y=se_y,
               centre=centre, bx_c=bx_c)
    regions.append(reg)
    genes.append((name, chrom, centre - 15_000, centre + 15_000, "+"))


for name, chrom, role, theta in PATHWAY:
    build(name, chrom, role, theta)
for k in range(N_GW):
    build(f"SYNG{k + 1:02d}", 1 + (k % 22), "gw", THETA_GW, is_gw=True)

# PLINK requires each chromosome to be contiguous in the .bim: sort regions by (chr, pos).
regions.sort(key=lambda r: (r["chr"], int(r["pos"][0])))

# --- flatten to SNP tables -------------------------------------------------
rows = []          # one dict per SNP, in region order (also the .bim/.bed order)
sid = 0
for r_i, reg in enumerate(regions):
    for i in range(L):
        sid += 1
        rows.append(dict(id=f"s{sid:06d}", region=r_i, i=i, chr=reg["chr"], pos=int(reg["pos"][i]),
                         a1=reg["a1"][i], a2=reg["a2"][i], eaf=float(reg["eaf"][i]),
                         bx=float(reg["bx"][i]), se_x=float(reg["se_x"][i]),
                         by=float(reg["by"][i]), se_y=float(reg["se_y"][i])))

# BCAA GWAS: original orientation
with gzip.open(DATA / "bcaa_gwas.tsv.gz", "wt", newline="\n") as f:
    f.write("SNP\tCHR\tBP\tA1\tA2\tFREQ_A1\tBETA\tSE\tP\tN\n")
    for r in rows:
        f.write(f"{r['id']}\t{r['chr']}\t{r['pos']}\t{r['a1']}\t{r['a2']}\t{r['eaf']:.4f}\t"
                f"{r['bx']:.5f}\t{r['se_x']:.5f}\t{pval(r['bx'] / r['se_x']):.4g}\t{N_BCAA}\n")

# T2D GWAS: different columns, shuffled rows, 15% allele-flipped
perm = rng.permutation(len(rows))
flip = rng.random(len(rows)) < 0.15
with gzip.open(DATA / "t2d_gwas.tsv.gz", "wt", newline="\n") as f:
    f.write("variant_id\tchromosome\tposition\teffect_allele\tother_allele\teffect_allele_frequency\t"
            "beta\tstandard_error\tp_value\tn_cases\tn_controls\n")
    for k in perm:
        r = rows[k]
        ea, oa, fr, b = r["a1"], r["a2"], r["eaf"], r["by"]
        if flip[k]:
            ea, oa, fr, b = oa, ea, 1 - fr, -b
        f.write(f"{r['id']}\t{r['chr']}\t{r['pos']}\t{ea}\t{oa}\t{fr:.4f}\t{b:.5f}\t{r['se_y']:.5f}\t"
                f"{pval(r['by'] / r['se_y']):.4g}\t{N_CASE}\t{N_CTRL}\n")

# gene annotation (fictional coordinates on the same synthetic build)
with open(DATA / "genes.tsv", "w", newline="\n") as f:
    f.write("gene\tchr\tstart\tend\tstrand\n")
    for g in genes:
        f.write("\t".join(map(str, g)) + "\n")

# PLINK 1 binary reference panel (500 individuals, same SNP ids, A1 = counted allele)
with open(DATA / "ld_ref.bim", "w", newline="\n") as f:
    for r in rows:
        f.write(f"{r['chr']}\t{r['id']}\t0\t{r['pos']}\t{r['a1']}\t{r['a2']}\n")
with open(DATA / "ld_ref.fam", "w", newline="\n") as f:
    for j in range(N_REF):
        f.write(f"REF{j:03d} REF{j:03d} 0 0 0 -9\n")
with open(DATA / "ld_ref.bed", "wb") as f:
    f.write(bytes([0x6C, 0x1B, 0x01]))
    code = {2: 0b00, 1: 0b10, 0: 0b11}        # count of A1 -> plink 2-bit code
    pad = (-N_REF) % 4
    for reg in regions:
        for i in range(L):
            g = reg["ref"][:, i]
            codes = np.array([code[int(x)] for x in g] + [0b01] * pad, dtype=np.uint8).reshape(-1, 4)
            f.write((codes[:, 0] | (codes[:, 1] << 2) | (codes[:, 2] << 4) | (codes[:, 3] << 6)).astype(np.uint8).tobytes())

# --- eval.json ---------------------------------------------------------------
TASK = """Your lab is investigating whether lowering plasma BCAA levels via inhibition of BCKDK kinase could reduce T2D risk.

In order to do so, you have obtained GWAS summary statistics on T2D risk, as well as circulating BCAA levels, all in European ancestry populations (data/bcaa_gwas.tsv.gz, data/t2d_gwas.tsv.gz, gene coordinates in data/genes.tsv, and a PLINK LD reference panel data/ld_ref.*).

As the degree of genetic correlation between circulating metabolites and heterogeneity between variants suggests a high degree of horizontal pleiotropy and difficulty in interpretation, you wish to avoid these issues while modeling the predicted change in T2D risk per 1 SD increase in levels of circulating BCAA levels at 3 best selected loci.

What is the highest change you observe and the lowest across the loci that you investigate (report estimates as log odds ratios per 1 SD increase in BCAA)? Does genetic evidence support the idea that BCKDK inhibition would be beneficial for reducing T2D risk (answer "YES" or "NO") in a list?

Return EXACTLY:

<EVAL_ANSWER>
{"loci_gene_names": [], "highest_estimate": "", "lowest_estimate": "", "bckdk_inhibition_yn": []}
</EVAL_ANSWER>"""

# ground truth is filled in by solvers.py after the reference estimator runs
truth = {"seed": SEED, "theta_true": {n: t for n, _, role, t in PATHWAY if role == "signal"},
         "signal_genes": [n for n, _, role, _ in PATHWAY if role == "signal"]}
(HERE / "truth.json").write_text(json.dumps(truth, indent=2), encoding="utf-8")
EVAL = {
    "id": "synthetic_cis_mr_bcaa",
    "replicates": "variantbench/cis_mr_bcaa",
    "synthetic": True,
    "task": TASK,
    "graders": "PLACEHOLDER - written by solvers.py --write-eval",
    "data_dir": "data",
}
(HERE / "eval.json").write_text(json.dumps(EVAL, indent=2), encoding="utf-8")
print("SNPs:", len(rows), "regions:", len(regions))
for reg in regions[:8]:
    print(f"  {reg['name']:7s} chr{reg['chr']:<2d} {reg['role']:6s} max|z_x|={np.abs(reg['bx'] / reg['se_x']).max():.1f}")
