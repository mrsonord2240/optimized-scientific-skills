"""Synthetic label-free proteomics differential-abundance task (i6_proteomics_de).

Everything here is FABRICATED from a known model.

Model (log2 scale):
  y[i,j] = mu_i + batch_i*[batch_j==2] + fc_i*[treated_j] + load_j + e_ij
  mu_i        ~ N(MU0, MU_SD)                       (log-normal-ish protein abundance)
  e_ij        ~ N(0, s_i^2), s_i^2 = s0(mu_i)^2 * D0 / chi2_D0   (noisier at low abundance)
  load_j      : per-sample loading offset, range ~1 log2 (up to ~2x), mean-imbalanced across
                the two conditions by chance (rejection-sampled, see below)
  batch_i     ~ N(0, BATCH_SD) per-protein shift for batch 2; 2 batches, 3+3 in each condition
  fc_i        : 150 proteins, |log2FC| U(0.9, 2.0), 50/50 up/down, positions random
Missingness (MNAR + a little MCAR): each cell is dropped with prob expit(-(y - X0)/S), i.e. the
probability of dropout rises steeply as the (sample-loaded) intensity nears the detection limit,
plus 2% random dropout. Missing cells are written as empty strings.
Table is MaxQuant-like: raw linear-scale intensities, Reverse / Potential contaminant flags.

Truth = the true-changed proteins that are *eligible* (>= MIN_VALID valid values in each condition);
that eligibility rule is stated in the task text, so it is an observable property of the data.

Loading offsets are drawn by rejection so the treated-minus-control mean offset is 0.30-0.45 log2 in
magnitude: a chance imbalance of the size a real 6-vs-6 experiment can show, which makes un-normalised
comparisons visibly biased. Fixed seed, so fully deterministic.

Run: python generate.py -> data/, truth.json, eval.json (graders written here; solvers.py checks them)
"""
import json
import os
from pathlib import Path

import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
SEED = int(os.environ.get("SEED", 1))
rng = np.random.default_rng(SEED)

N_REAL, N_REV, N_CON = 1480, 12, 8
N_DE = 150
MU0, MU_SD = 23.3, 2.4
X0, S_DROP = 21.4, 0.55
P_MCAR = 0.02
BATCH_SD = 0.6
D0 = 8.0
MIN_VALID = 3
THRESH = {"precision_at_k": 0.85, "recall_at_k": 0.80}   # calibrated, see NOTES.md

SAMPLES = ["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11", "S12"]
COND = ["treated", "control", "control", "treated", "treated", "control",
        "control", "treated", "control", "treated", "treated", "control"]
BATCH = [1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 2]   # 3 treated + 3 control in each batch


def build():
    cond = np.array(COND)
    batch = np.array(BATCH)
    treated = (cond == "treated").astype(float)
    b2 = (batch == 2).astype(float)
    n = len(SAMPLES)

    # loading offsets: uniform over ~1 log2 (2x); centred overall
    while True:
        load = rng.uniform(-0.5, 0.5, n)
        load -= load.mean()
        if 0.30 <= abs(load[treated == 1].mean() - load[treated == 0].mean()) <= 0.45:
            break

    n_all = N_REAL + N_REV + N_CON
    mu = rng.normal(MU0, MU_SD, n_all)
    # sample-specific noise sd, larger at low abundance, moderated inverse-chi2 spread
    s0 = 0.20 + 0.45 * expit(-(mu - 21.0) / 1.3)
    s2 = (s0 ** 2) * D0 / rng.chisquare(D0, n_all)
    sd = np.sqrt(s2)
    bshift = rng.normal(0, BATCH_SD, n_all)

    fc = np.zeros(n_all)
    de_idx = rng.choice(N_REAL, N_DE, replace=False)
    mag = rng.uniform(0.9, 2.0, N_DE)
    sign = np.where(np.arange(N_DE) % 2 == 0, 1.0, -1.0)
    rng.shuffle(sign)
    fc[de_idx] = mag * sign
    # a few contaminants/reverse hits with large condition differences (keratin-like carry-over)
    con_idx = np.arange(N_REAL + N_REV, n_all)
    fc[con_idx[:4]] = np.array([1.6, 1.3, -1.4, 1.8])

    y = (mu[:, None] + bshift[:, None] * b2[None, :] + fc[:, None] * treated[None, :]
         + load[None, :] + rng.normal(0, 1, (n_all, n)) * sd[:, None])

    p_drop = expit(-(y - X0) / S_DROP)
    p_drop = 1 - (1 - p_drop) * (1 - P_MCAR)
    miss = rng.random((n_all, n)) < p_drop
    return dict(y=y, miss=miss, mu=mu, fc=fc, de_idx=de_idx, load=load, cond=cond, batch=batch)


def main():
    w = build()
    y, miss, fc = w["y"], w["miss"], w["fc"]
    n_all = y.shape[0]
    raw = np.where(miss, np.nan, np.round(2.0 ** y))

    # identifiers
    letters = np.array(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
    ids = set()
    pid, gene = [], []
    while len(pid) < n_all:
        a = "".join(rng.choice(list("OPQ"), 1)) + str(rng.integers(0, 10)) + "".join(rng.choice(letters, 3)) + str(rng.integers(0, 10))
        if a in ids:
            continue
        ids.add(a)
        pid.append(a)
        gene.append("".join(rng.choice(letters, 3)) + str(rng.integers(1, 40)) + (rng.choice(["", "A", "B"]) if rng.random() < .3 else ""))
    # ensure unique gene names
    seen = {}
    for i, g in enumerate(gene):
        if g in seen:
            g = g + "L"
        while g in seen:
            g = g + "X"
        seen[g] = 1
        gene[i] = g
    reverse = np.zeros(n_all, bool)
    reverse[N_REAL:N_REAL + N_REV] = True
    contam = np.zeros(n_all, bool)
    contam[N_REAL + N_REV:] = True
    pid = [("REV__" + p) if reverse[i] else ("CON__" + p) if contam[i] else p for i, p in enumerate(pid)]
    gene = ["" if reverse[i] else g for i, g in enumerate(gene)]

    # peptide counts, correlated with abundance
    mu = w["mu"]
    pept = np.clip(np.round(1 + np.maximum(mu - 20, 0) * 1.6 + rng.normal(0, 2, n_all)), 1, 60).astype(int)

    cols = ["Protein IDs", "Gene names", "Peptides", "Reverse", "Potential contaminant"] + \
           [f"Intensity {s}" for s in SAMPLES]
    lines = ["\t".join(cols)]
    order = rng.permutation(n_all)          # shuffle rows so decoys/contaminants are not clustered
    for i in order:
        cells = [pid[i], gene[i], str(pept[i]), "+" if reverse[i] else "", "+" if contam[i] else ""]
        cells += ["" if np.isnan(v) else str(int(v)) for v in raw[i]]
        lines.append("\t".join(cells))
    (DATA / "protein_intensities.tsv").write_bytes(("\n".join(lines) + "\n").encode("utf-8"))

    sl = ["sample\tcondition\tbatch"] + [f"{s}\t{c}\t{b}" for s, c, b in zip(SAMPLES, COND, BATCH)]
    (DATA / "samples.tsv").write_bytes(("\n".join(sl) + "\n").encode("utf-8"))

    # truth: true-changed real proteins that are eligible
    cond = w["cond"]
    valid = ~miss
    ok_t = valid[:, cond == "treated"].sum(1) >= MIN_VALID
    ok_c = valid[:, cond == "control"].sum(1) >= MIN_VALID
    eligible = ok_t & ok_c
    de_all = [pid[i] for i in w["de_idx"]]
    truth_ids = sorted(pid[i] for i in w["de_idx"] if eligible[i])
    truth = {"seed": SEED, "n_de_designed": int(N_DE), "n_de_eligible": len(truth_ids),
             "de_all": {pid[i]: float(fc[i]) for i in w["de_idx"]},
             "de_eligible": truth_ids,
             "loading_log2": dict(zip(SAMPLES, [float(x) for x in w["load"]])),
             "mean_missing_fraction": float(miss[:N_REAL].mean()),
             "note": "truth = designed-changed proteins with >= %d valid values in each condition" % MIN_VALID}
    (HERE / "truth.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")

    # eval.json: thresholds calibrated in solvers.py / NOTES.md
    ev_path = HERE / "eval.json"
    graders = [{"type": "marker_gene_precision_recall",
                "config": {"canonical_markers": truth_ids, "answer_field": "differentially_abundant_proteins",
                           "scoring": {"pass_thresholds": THRESH}}}]
    ev = {"id": "i6_proteomics_de", "replicates": None, "synthetic": True,
          "task": TASK, "graders": graders, "data_dir": "data",
          "relevant_skills": RELEVANT,
          "design_source": "Label-free proteomics missing values are largely missing-not-at-random (low-abundance "
                           "proteins drop out) and per-sample loading differs; naive min/zero imputation plus t-test "
                           "manufactures false positives (Lazar et al. 2016, J Proteome Res 15:1116)."}
    ev_path.write_text(json.dumps(ev, indent=2), encoding="utf-8")
    print(f"wrote data/, truth.json, eval.json | seed={SEED} designed DE={N_DE} eligible DE={len(truth_ids)} "
          f"missing={truth['mean_missing_fraction']:.3f} load(treated-control)="
          f"{w['load'][cond=='treated'].mean()-w['load'][cond=='control'].mean():+.3f}")


TASK = (
    "A collaborator ran a label-free quantitative proteomics experiment comparing 6 treated and 6 control samples "
    "(processed in two batches, each containing 3 treated and 3 control samples). The protein-level output of the search "
    "software is in data/protein_intensities.tsv: one row per protein group, raw (linear-scale) intensities in the columns "
    "named `Intensity <sample>`, empty cells where a protein was not quantified in a sample, plus the columns "
    "`Reverse` and `Potential contaminant` (a '+' marks a flagged row). Sample metadata (condition, batch) are in "
    "data/samples.tsv; batch is a known nuisance variable.\n\n"
    "Identify the proteins that differ in abundance between treated and control, controlling the false discovery rate at 5% "
    "with the Benjamini-Hochberg procedure across all proteins you test. Do not report reverse-database hits or potential "
    "contaminants. Only proteins with at least 3 valid (non-missing) intensity values in each condition are eligible to be "
    "reported. Direction of change does not matter; report a protein if it is significantly different in either direction.\n\n"
    "Return EXACTLY:\n\n<EVAL_ANSWER>\n{\"differentially_abundant_proteins\": [\"<Protein IDs value>\", ...]}\n</EVAL_ANSWER>\n\n"
    "where each list item is the value in the `Protein IDs` column."
)
RELEVANT = ["bio-proteomics-differential-abundance", "bio-proteomics-quantification", "bio-proteomics-data-import",
            "bio-proteomics-proteomics-qc", "bio-experimental-design-multiple-testing",
            "bio-experimental-design-batch-design"]

if __name__ == "__main__":
    main()
