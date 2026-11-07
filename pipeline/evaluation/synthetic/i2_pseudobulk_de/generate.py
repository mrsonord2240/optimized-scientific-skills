"""Synthetic eval task i2_pseudobulk_de.

Everything here is FABRICATED. Multi-donor scRNA-seq of ONE annotated cell type, control vs treated,
each donor sampled in both conditions with very unequal cell counts. Counts are drawn from a
negative-binomial model with
  * donor-specific baseline shifts (large on a subset of genes, small elsewhere),
  * a donor x condition sample-level random effect (biological noise between samples),
  * cell-level gamma-Poisson noise and library-size variation.
A known set of genes carries a true treatment effect |log2FC| >= 1.4.

Run: python generate.py -> data/, truth.json, eval.json   (deterministic, seed 1)
"""
import gzip
import io
import json
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
SEED = int(os.environ.get("SEED", 1))

G = 4000
N_DE = 300                 # true DE genes, |log2FC| in [1.4, 2.6]
N_MODEST = 150             # real but small effects, |log2FC| in [0.3, 0.7] (not in truth)
V_FRAC = 0.20              # fraction of genes with large donor effects
SIG_DONOR_V, SIG_DONOR_BASE = 1.5, 0.60     # natural-log SD of donor effect
SIG_SAMPLE = 0.20          # natural-log SD of donor x condition effect (all genes)
SIZE_MEAN, SIZE_SD = 6.0, 0.4   # NB size (1/dispersion), lognormal across genes
DEPTH = 2200.0             # mean UMI per cell
MIN_CELLS = 5              # donors with fewer cells in either condition are excluded (task text)
MIN_TOTAL_COUNT = 20       # gene filter (task text)

# (control cells, treated cells) per donor; donors shuffled to opaque ids below
CELLS = [(300, 40), (340, 55), (260, 25), (110, 280), (80, 340), (60, 240), (140, 180), (45, 12),
         (150, 2), (220, 3), (4, 300), (180, 1)]


def build():
    rng = np.random.default_rng(SEED)
    nd = len(CELLS)

    # ---- gene-level model --------------------------------------------------------------
    logp = rng.normal(0.0, 1.6, G)
    p = np.exp(logp)
    p /= p.sum()
    mean_cell = DEPTH * p
    size_g = np.exp(rng.normal(np.log(SIZE_MEAN), SIZE_SD, G))

    idx_ok = np.where(mean_cell >= 0.5)[0]
    perm = rng.permutation(idx_ok)
    de = np.sort(perm[:N_DE])
    modest = np.sort(perm[N_DE:N_DE + N_MODEST])
    lfc = np.zeros(G)                               # true log2FC treated vs control
    # Signs are chosen greedily so that up- and down-regulated expression mass balances (no compositional
    # shift of the null genes under library-size normalisation): the truth is then normaliser-independent.
    mags = {int(i): m for i, m in zip(de, rng.uniform(1.4, 2.6, N_DE))}
    mags.update({int(i): m for i, m in zip(modest, rng.uniform(0.3, 0.7, N_MODEST))})
    S = 0.0
    for i in sorted(mags, key=lambda i: -p[i]):
        up, dn = p[i] * (2 ** mags[i] - 1), p[i] * (2 ** -mags[i] - 1)
        if abs(S + up) <= abs(S + dn):
            lfc[i], S = mags[i], S + up
        else:
            lfc[i], S = -mags[i], S + dn
    print("composition drift of treated fractions before renorm: log2 R =", round(float(np.log2(1 + S)), 4))

    is_v = rng.random(G) < V_FRAC
    sig_d = np.where(is_v, SIG_DONOR_V, SIG_DONOR_BASE)
    U = rng.normal(0.0, 1.0, (nd, G)) * sig_d                       # donor effect (ln)
    E = rng.normal(0.0, SIG_SAMPLE, (nd, 2, G))                     # donor x condition (ln)
    donor_depth = np.exp(rng.normal(0.0, 0.2, (nd, 2)))

    # ---- cells --------------------------------------------------------------------------
    rows_d, rows_c = [], []
    for d, (nc, nt) in enumerate(CELLS):
        rows_d += [d] * (nc + nt)
        rows_c += [0] * nc + [1] * nt
    rows_d, rows_c = np.array(rows_d), np.array(rows_c)
    n = len(rows_d)
    lib = DEPTH * donor_depth[rows_d, rows_c] * np.exp(rng.normal(0.0, 0.35, n) - 0.35 ** 2 / 2)
    counts = np.zeros((n, G), dtype=np.int32)
    for d in range(nd):
        for c in (0, 1):
            sel = np.where((rows_d == d) & (rows_c == c))[0]
            if len(sel) == 0:
                continue
            lg = logp - np.log(np.exp(logp).sum()) + U[d] + E[d, c] + c * lfc * np.log(2)
            frac = np.exp(lg)
            frac /= frac.sum()                       # compositional: relative abundance within cell
            mu = lib[sel, None] * frac[None, :]
            lam = mu * rng.gamma(size_g, 1.0 / size_g, (len(sel), G))
            counts[sel] = rng.poisson(lam)

    # shuffle cell order; opaque donor ids and barcodes
    order = rng.permutation(n)
    counts, rows_d, rows_c = counts[order], rows_d[order], rows_c[order]
    donor_ids = rng.permutation(np.arange(101, 101 + nd))
    donor_name = np.array([f"D{donor_ids[d]}" for d in range(nd)])
    bases = np.array(list("ACGT"))
    bc = set()
    while len(bc) < n:
        bc.add("".join(rng.choice(bases, 16)))
    bc = sorted(bc)
    rng.shuffle(bc)
    barcodes = [b + "-1" for b in bc]
    return dict(counts=counts, donor=donor_name[rows_d], cond=np.where(rows_c == 1, "treated", "control"),
                barcodes=barcodes, de=de, modest=modest, is_v=is_v, lfc=lfc, U=U, mean_cell=mean_cell)


def gz_write(path, text):
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0, compresslevel=6) as f:
        f.write(text.encode("utf-8"))
    Path(path).write_bytes(buf.getvalue())


def main():
    w = build()
    counts = w["counts"]
    n, g = counts.shape
    ids = [f"ENSGSYN{i + 1:07d}" for i in range(g)]
    syms = [f"SYNG{i + 1:04d}" for i in range(g)]
    gi, ci = np.nonzero(counts.T)                    # column-major (cell by cell)
    v = counts.T[gi, ci]
    lines = ["%%MatrixMarket matrix coordinate integer general", f"{g} {n} {len(v)}"]
    lines += [f"{a + 1} {b + 1} {c}" for a, b, c in zip(gi, ci, v)]
    gz_write(DATA / "matrix.mtx.gz", "\n".join(lines) + "\n")
    gz_write(DATA / "features.tsv.gz", "".join(f"{i}\t{s}\tGene Expression\n" for i, s in zip(ids, syms)))
    gz_write(DATA / "barcodes.tsv.gz", "".join(b + "\n" for b in w["barcodes"]))
    with open(DATA / "cells.tsv", "w", newline="\n") as f:
        f.write("barcode\tdonor\tcondition\tcell_type\n")
        for b, d, c in zip(w["barcodes"], w["donor"], w["cond"]):
            f.write(f"{b}\t{d}\t{c}\tCD4 T cell\n")

    truth_genes = [syms[i] for i in w["de"]]
    truth = {"de_genes": truth_genes,
             "true_log2fc": {syms[i]: round(float(w["lfc"][i]), 4) for i in np.where(w["lfc"] != 0)[0]},
             "modest_effect_genes": [syms[i] for i in w["modest"]],
             "large_donor_effect_genes": [syms[i] for i in np.where(w["is_v"])[0]]}
    (HERE / "truth.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")

    task = (
        "data/ holds a single-cell RNA-seq experiment of one already-annotated cell type (CD4 T cells) collected from "
        "12 donors: raw integer UMI counts in 10x format (data/matrix.mtx.gz with genes as rows and cells as columns, "
        "data/features.tsv.gz, data/barcodes.tsv.gz) and data/cells.tsv giving each cell's donor and condition "
        "(control or treated). Each donor was sampled in both conditions, but the number of cells per donor and "
        "condition varies widely.\n\n"
        "Find the genes that respond to treatment (treated vs control). Report every gene that meets ALL of:\n"
        "  1. Adjusted p-value < 0.05 (Benjamini-Hochberg across all tested genes).\n"
        "  2. |log2 fold change| > 1, where log2 fold change is the log2 of the ratio of mean library-size-normalised "
        "expression (normalised-count scale, not a mean of log values) in treated over control.\n"
        "Rules for the analysis: (a) drop entirely (in both conditions) any donor that has fewer than 5 cells in "
        "control or fewer than 5 cells in treated; (b) test only genes with a total raw count of at least "
        f"{MIN_TOTAL_COUNT} summed over all cells of the retained donors.\n\n"
        "Return the gene symbols (column 2 of data/features.tsv.gz) in this exact format:\n\n"
        "<EVAL_ANSWER>\n{\"de_genes\": [\"<symbol>\", \"<symbol>\", ...]}\n</EVAL_ANSWER>"
    )
    ev = {"id": "i2_pseudobulk_de", "replicates": None, "synthetic": True, "task": task,
          "graders": [{"type": "marker_gene_precision_recall", "config": {
              "canonical_markers": truth_genes, "answer_field": "de_genes",
              "scoring": {"pass_thresholds": {"precision_at_k": PREC, "recall_at_k": REC}}}}],
          "data_dir": "data",
          "relevant_skills": RELEVANT_SKILLS,
          "design_source": "Treating cells as independent replicates inflates significance in multi-donor scRNA-seq DE; "
                           "donor is the unit of replication (Squair et al. 2021 Nat Commun; Murphy & Skene 2022 Nat Commun)."}
    (HERE / "eval.json").write_text(json.dumps(ev, indent=2), encoding="utf-8")
    print(f"cells={n} genes={g} nnz={len(v)} de={len(truth_genes)}")


PREC, REC = 0.65, 0.70
RELEVANT_SKILLS = ["bio-differential-expression-deseq2-basics", "bio-single-cell-preprocessing", "bio-single-cell-markers-annotation",
                   "bio-single-cell-data-io", "bio-workflows-scrnaseq-pipeline", "bio-experimental-design-randomization-blocking"]

if __name__ == "__main__":
    main()
