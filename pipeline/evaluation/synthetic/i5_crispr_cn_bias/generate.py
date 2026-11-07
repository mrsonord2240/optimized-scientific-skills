"""Synthetic pooled CRISPR-Cas9 knockout screen in one cancer cell line.

Everything here is FABRICATED (gene symbols invented, SYN####; coordinates, copy
number, counts all simulated).  Ground truth is known by construction.

Model (per guide, log2 fold change late vs plasmid):
    L = e_guide * (E_gene + C_cn)
  E_gene : true gene fitness effect  (strong essential, weak, or neutral)
  C_cn   : copy-number-dependent cutting effect  (gain: -1.2*log2(CN/2);
           loss: +0.3*log2(2/CN)), independent of the gene's function
  e      : guide efficiency ~ Beta(4, 1.5)
Design source: Aguirre et al. 2016 Cancer Discov; Meyers et al. 2017 Nat Genet
(copy-number-dependent Cas9 cutting toxicity confounds knockout screens).

Run:  python generate.py  -> data/*.tsv, truth.json, eval.json
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
rng = np.random.default_rng(20260919)

# ---- knobs -----------------------------------------------------------------
N_GENES = 600
GUIDES_PER_GENE = 4
N_NTC = 100
N_ESS = 40           # true strong dependencies (2 of them sit inside the amplicon)
N_WEAK = 20          # weakly depleting
SLOPE_GAIN = 1.2     # log2FC per log2(CN/2) of gain
SHIFT_LOSS = 0.3     # positive shift per log2(2/CN) of loss
PRECISION_MIN = 0.75
RECALL_MIN = 0.75

CHROMS = ["chr%d" % i for i in range(1, 13)]
LEN_MB = [140, 125, 110, 100, 95, 90, 85, 80, 70, 65, 60, 55]
CHR_LEN = {c: l * 1_000_000 for c, l in zip(CHROMS, LEN_MB)}
AMP_CHR, AMP_WIN = "chr7", (40_000_000, 50_000_000)
AMP_GROUPS = [(5, 5.2), (6, 7.0), (9, 10.5), (6, 8.6), (6, 5.6)]   # (n genes, copy number)
N_AMP = sum(n for n, _ in AMP_GROUPS)                               # 32 = 30 non-ess + 2 ess

# ---- genes and coordinates -------------------------------------------------
tot = sum(LEN_MB)
n_per = {c: max(30, round((N_GENES - 55) * l / (tot - 90))) for c, l in zip(CHROMS, LEN_MB)}
n_per[AMP_CHR] = 55
n_per["chr1"] += N_GENES - sum(n_per.values())
assert sum(n_per.values()) == N_GENES


def place(chrom, n, lo, hi):
    starts = np.sort(rng.integers(lo, hi, n))
    out = []
    for i, s in enumerate(starts):
        gap = (starts[i + 1] - s) if i + 1 < n else (hi - s)
        ln = int(np.clip(rng.lognormal(10.6, 0.6), 5_000, 120_000))
        out.append((chrom, int(s), int(s + min(ln, max(2_000, int(gap * 0.5))))))
    return out


genes = []          # dicts in chromosome order, sorted by start within chromosome
n_out = n_per[AMP_CHR] - N_AMP
n_a = n_out // 2    # number of chr7 genes upstream of the amplicon window
for c in CHROMS:
    if c == AMP_CHR:
        pos = (place(c, n_a, 1_000_000, AMP_WIN[0] - 2_000_000)
               + place(c, N_AMP, AMP_WIN[0], AMP_WIN[1])
               + place(c, n_out - n_a, AMP_WIN[1] + 2_000_000, CHR_LEN[c] - 1_000_000))
    else:
        pos = place(c, n_per[c], 1_000_000, CHR_LEN[c] - 1_000_000)
    for i, (ch, s, e) in enumerate(pos):
        genes.append(dict(chr=ch, start=s, end=e, idx=i))
for g in genes:
    g["amp"] = (g["chr"] == AMP_CHR and n_a <= g["idx"] < n_a + N_AMP)

names = ["SYN%04d" % k for k in rng.choice(np.arange(1, 10000), N_GENES, replace=False)]
for g, nm in zip(genes, names):
    g["gene"] = nm

# ---- copy-number blocks -> segments ----------------------------------------
blocks_def = {}
n3 = n_per["chr3"]
blocks_def["chr3"] = [(0, int(0.55 * n3) - 1, 3.0)]       # arm-level gain
n5 = n_per["chr5"]
blocks_def["chr5"] = [(0, int(0.45 * n5) - 1, 1.0)]       # arm-level loss
blocks_def["chr10"] = [(10, 17, 4.0)]                      # small focal gain
blocks_def["chr12"] = [(3, 8, 1.0)]                        # small focal loss
amp_blocks, a = [], n_a
for n, cn in AMP_GROUPS:
    amp_blocks.append((a, a + n - 1, cn))
    a += n
blocks_def[AMP_CHR] = amp_blocks

seg_rows = []
gene_cn = {}
for c in CHROMS:
    cg = sorted([g for g in genes if g["chr"] == c], key=lambda g: g["start"])
    n = len(cg)
    lab = [2.0] * n
    bid = [-1] * n
    for k, (i0, i1, cn) in enumerate(blocks_def.get(c, [])):
        for i in range(i0, i1 + 1):
            lab[i], bid[i] = cn, k
    # runs of identical (block id, cn); long diploid runs are split in two
    runs, i = [], 0
    while i < n:
        j = i
        while j + 1 < n and bid[j + 1] == bid[i] and lab[j + 1] == lab[i]:
            j += 1
        if bid[i] == -1 and (j - i + 1) > 12:
            m = i + int(rng.integers(4, j - i - 3))
            runs += [(i, m, 2.0), (m + 1, j, 2.0)]
        else:
            runs.append((i, j, lab[i]))
        i = j + 1
    for r, (i0, i1, cn) in enumerate(runs):
        lo = 0 if r == 0 else (cg[i0 - 1]["end"] + cg[i0]["start"]) // 2
        hi = CHR_LEN[c] if r == len(runs) - 1 else (cg[i1]["end"] + cg[i1 + 1]["start"]) // 2
        rep = round(cn + (rng.normal(0, 0.025) if cn == 2.0 else rng.normal(0, 0.02) * cn), 3)
        seg_rows.append((c, lo, hi, rep))
        for g in cg[i0:i1 + 1]:
            gene_cn[g["gene"]] = rep
for g in genes:
    g["cn"] = gene_cn[g["gene"]]

# ---- roles -------------------------------------------------------------------
amp_genes = [g for g in genes if g["amp"]]
other = [g for g in genes if not g["amp"]]
ess = [amp_genes[i] for i in rng.choice(N_AMP, 2, replace=False)]
ess += [other[i] for i in rng.choice(len(other), N_ESS - 2, replace=False)]
ess_set = {g["gene"] for g in ess}
rest = [g for g in other if g["gene"] not in ess_set]
weak = [rest[i] for i in rng.choice(len(rest), N_WEAK, replace=False)]
weak_set = {g["gene"] for g in weak}
for g in genes:
    if g["gene"] in ess_set:
        g["E"] = -float(rng.uniform(2.4, 5.0))
    elif g["gene"] in weak_set:
        g["E"] = -float(rng.uniform(0.3, 0.6))
    else:
        g["E"] = 0.0
    x = np.log2(g["cn"] / 2.0)
    g["C"] = -SLOPE_GAIN * x if x > 0 else SHIFT_LOSS * (-x)

# ---- guides ------------------------------------------------------------------
BASES = np.array(list("ACGT"))
guides = []   # (guide_id, gene, seq, L)
for g in sorted(genes, key=lambda g: g["gene"]):
    for k in range(GUIDES_PER_GENE):
        e = float(rng.beta(4, 1.5))
        L = max(-7.0, e * (g["E"] + g["C"]))
        seq = "".join(BASES[rng.integers(0, 4, 20)])
        guides.append((f"{g['gene']}_sg{k + 1}", g["gene"], seq, L))
for k in range(N_NTC):
    seq = "".join(BASES[rng.integers(0, 4, 20)])
    guides.append((f"NTC_{k + 1:04d}", "NTC", seq, 0.0))
G = len(guides)
L = np.array([x[3] for x in guides])

# ---- counts ------------------------------------------------------------------
p = np.exp(rng.normal(0, 0.5, G))
p /= p.sum()


def nb(mu, size):
    lam = rng.gamma(size, mu / size)
    return rng.poisson(lam)


plasmid = nb(600 * G * p, 30)


def late(depth):
    q = p * 2 ** (L + rng.normal(0, 0.2, G))
    q /= q.sum()
    return nb(depth * G * q, 40)


rep1, rep2 = late(500), late(450)


# ---- write -------------------------------------------------------------------
def w(path, header, rows):
    with open(path, "w", newline="\n") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(v) for v in r) + "\n")


w(DATA / "counts.tsv", ["guide_id", "plasmid", "late_rep1", "late_rep2"],
  [(guides[i][0], plasmid[i], rep1[i], rep2[i]) for i in range(G)])
w(DATA / "library.tsv", ["guide_id", "gene", "sequence"], [(x[0], x[1], x[2]) for x in guides])
gs = sorted(genes, key=lambda g: (CHROMS.index(g["chr"]), g["start"]))
w(DATA / "gene_coordinates.tsv", ["gene", "chromosome", "start", "end"],
  [(g["gene"], g["chr"], g["start"], g["end"]) for g in gs])
w(DATA / "copy_number_segments.tsv", ["chromosome", "start", "end", "copy_number"], seg_rows)

truth = {
    "essential_genes": sorted(ess_set),
    "essential_in_amplicon": sorted(g["gene"] for g in ess if g["amp"]),
    "weak_genes": sorted(weak_set),
    "amplicon_non_essential": sorted(g["gene"] for g in amp_genes if g["gene"] not in ess_set),
    "gene_cn": {g["gene"]: g["cn"] for g in genes},
    "gene_E": {g["gene"]: round(g["E"], 3) for g in genes},
}
(HERE / "truth.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")

TASK = (
    "We ran a pooled CRISPR-Cas9 knockout screen in a single cancer cell line (library of 4 sgRNAs per gene plus "
    "non-targeting control guides, labelled NTC). Cells were sampled at the plasmid stage and again at a late "
    "timepoint in two independent replicates. Files:\n"
    "- data/counts.tsv: raw sgRNA read counts (rows = guides; columns: plasmid, late_rep1, late_rep2)\n"
    "- data/library.tsv: guide_id -> gene (and spacer sequence)\n"
    "- data/gene_coordinates.tsv: gene, chromosome, start, end (genome coordinates, one row per targeted gene)\n"
    "- data/copy_number_segments.tsv: chromosome, start, end, copy_number (segment-level total copy number of this "
    "cell line)\n\n"
    "Please identify the genes that are genuine essential dependencies of this cell line. Convention: a gene counts "
    "as an essential dependency if knocking it out substantially reduces the line's fitness (a gene-level log2 fold "
    "change of guide abundance, late versus plasmid, of about -1 or stronger, i.e. at least ~2-fold depletion, "
    "statistically supported at 5% FDR) for reasons intrinsic to loss of that gene, rather than as an artefact of "
    "the experiment. Do not report genes with only modest or marginal depletion.\n\n"
    "Return the gene symbols exactly as they appear in data/library.tsv, in this format:\n\n"
    "<EVAL_ANSWER>\n{\"essential_genes\": [\"<gene>\", \"<gene>\", ...]}\n</EVAL_ANSWER>"
)

ev = {
    "id": "i5_crispr_cn_bias", "replicates": None, "synthetic": True, "task": TASK,
    "graders": [{"type": "marker_gene_precision_recall", "config": {
        "canonical_markers": sorted(ess_set), "answer_field": "essential_genes",
        "scoring": {"pass_thresholds": {"precision_at_k": PRECISION_MIN, "recall_at_k": RECALL_MIN}}}}],
    "data_dir": "data",
    "relevant_skills": ["bio-crispr-screens-copy-number-correction", "bio-crispr-screens-mageck-analysis",
                        "bio-crispr-screens-hit-calling", "bio-crispr-screens-bagel-essentiality",
                        "bio-crispr-screens-screen-qc", "bio-workflows-crispr-screen-pipeline"],
    "design_source": "Copy-number-dependent Cas9 cutting toxicity makes non-essential genes in amplicons look like "
                     "dependencies in CRISPR KO screens (Aguirre et al. 2016 Cancer Discov; Meyers et al. 2017 Nat Genet).",
}
(HERE / "eval.json").write_text(json.dumps(ev, indent=2), encoding="utf-8")
print("wrote", G, "guides;", len(genes), "genes;", len(seg_rows), "CN segments;", len(ess_set), "essentials")
