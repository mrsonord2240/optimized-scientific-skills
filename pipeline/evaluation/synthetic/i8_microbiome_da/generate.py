"""Generate the synthetic microbiome differential-abundance task.

World model (truth by construction)
-----------------------------------
ABSOLUTE abundance a[i,j] of taxon i in sample j (e.g. cells per gram):
    log a = log base_i + logFC_i * case_j + sexeff_i * female_j + biological noise
Only 5 taxa bloom in cases and 3 taxa decline; every other taxon has the SAME
absolute abundance distribution in cases and controls. Sequencing sees only the
composition p = a / sum(a): the blooms inflate the denominator, so all unchanged
taxa fall in relative abundance in cases.
Reads: counts ~ Multinomial(depth_j, p_j * Gamma(k)-noise)  (extra-Poisson dispersion).

Truth set = taxa with a changed absolute abundance (5 up + 3 down), nothing else.

Run:  python generate.py            (writes data/, truth.json, eval.json)
Then: python solvers.py            (checks that the solvers pass/fail as designed)
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
D = HERE / "data"

SEED = 20260919
T, N_CTRL, N_CASE = 200, 20, 20
N = N_CTRL + N_CASE
N_ALLZERO, N_NEARZERO = 6, 9
UP_FC = [10.0, 8.0, 7.0, 6.0, 5.0]        # absolute fold change in cases
DOWN_FC = [0.25, 0.30, 0.35]
UP_SHARE = (0.10, 0.20)                   # combined baseline share of the 5 bloomers
N_SEX = 10                                 # taxa affected by sex (nuisance, NOT in truth)
PRECISION_MIN, RECALL_MIN = 0.40, 0.75    # grader thresholds (see NOTES.md for calibration)
GAMMA_K = 25.0                             # extra-Poisson dispersion of the sequencing step

LINEAGES = [  # (phylum, class, order, family, genus) - fabricated assignment, real names
    ("Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides"),
    ("Bacteroidota", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella"),
    ("Bacteroidota", "Bacteroidia", "Bacteroidales", "Rikenellaceae", "Alistipes"),
    ("Bacteroidota", "Bacteroidia", "Bacteroidales", "Tannerellaceae", "Parabacteroides"),
    ("Bacteroidota", "Bacteroidia", "Bacteroidales", "Barnesiellaceae", "Barnesiella"),
    ("Bacillota", "Clostridia", "Lachnospirales", "Lachnospiraceae", "Blautia"),
    ("Bacillota", "Clostridia", "Lachnospirales", "Lachnospiraceae", "Roseburia"),
    ("Bacillota", "Clostridia", "Lachnospirales", "Lachnospiraceae", "Coprococcus"),
    ("Bacillota", "Clostridia", "Lachnospirales", "Lachnospiraceae", "Dorea"),
    ("Bacillota", "Clostridia", "Lachnospirales", "Lachnospiraceae", "Anaerostipes"),
    ("Bacillota", "Clostridia", "Lachnospirales", "Lachnospiraceae", "Lachnoclostridium"),
    ("Bacillota", "Clostridia", "Oscillospirales", "Ruminococcaceae", "Faecalibacterium"),
    ("Bacillota", "Clostridia", "Oscillospirales", "Ruminococcaceae", "Ruminococcus"),
    ("Bacillota", "Clostridia", "Oscillospirales", "Ruminococcaceae", "Subdoligranulum"),
    ("Bacillota", "Clostridia", "Oscillospirales", "Oscillospiraceae", "Oscillibacter"),
    ("Bacillota", "Clostridia", "Eubacteriales", "Eubacteriaceae", "Eubacterium"),
    ("Bacillota", "Clostridia", "Peptostreptococcales", "Peptostreptococcaceae", "Romboutsia"),
    ("Bacillota", "Clostridia", "Clostridiales", "Clostridiaceae", "Clostridium"),
    ("Bacillota", "Negativicutes", "Veillonellales", "Veillonellaceae", "Veillonella"),
    ("Bacillota", "Negativicutes", "Acidaminococcales", "Acidaminococcaceae", "Phascolarctobacterium"),
    ("Bacillota", "Negativicutes", "Selenomonadales", "Selenomonadaceae", "Megamonas"),
    ("Bacillota", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"),
    ("Bacillota", "Bacilli", "Lactobacillales", "Streptococcaceae", "Streptococcus"),
    ("Bacillota", "Bacilli", "Lactobacillales", "Enterococcaceae", "Enterococcus"),
    ("Bacillota", "Erysipelotrichia", "Erysipelotrichales", "Erysipelotrichaceae", "Holdemanella"),
    ("Bacillota", "Erysipelotrichia", "Erysipelotrichales", "Erysipelotrichaceae", "Turicibacter"),
    ("Actinomycetota", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Bifidobacterium"),
    ("Actinomycetota", "Coriobacteriia", "Coriobacteriales", "Coriobacteriaceae", "Collinsella"),
    ("Actinomycetota", "Coriobacteriia", "Eggerthellales", "Eggerthellaceae", "Eggerthella"),
    ("Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia-Shigella"),
    ("Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Klebsiella"),
    ("Pseudomonadota", "Deltaproteobacteria", "Desulfovibrionales", "Desulfovibrionaceae", "Bilophila"),
    ("Pseudomonadota", "Betaproteobacteria", "Burkholderiales", "Sutterellaceae", "Sutterella"),
    ("Verrucomicrobiota", "Verrucomicrobiae", "Verrucomicrobiales", "Akkermansiaceae", "Akkermansia"),
    ("Fusobacteriota", "Fusobacteriia", "Fusobacteriales", "Fusobacteriaceae", "Fusobacterium"),
    ("Bacillota", "Clostridia", "Christensenellales", "Christensenellaceae", "Christensenellaceae_R-7_group"),
]


def build(seed=SEED):
    rng = np.random.default_rng(seed)

    # ---- sample design -----------------------------------------------------
    case = np.r_[np.zeros(N_CTRL, int), np.ones(N_CASE, int)]
    # sex: mild imbalance between groups (14F/6M in cases, 7F/13M in controls)
    female = np.r_[rng.permutation([1] * 7 + [0] * 13), rng.permutation([1] * 14 + [0] * 6)]
    order = rng.permutation(N)                         # shuffle sample order
    case, female = case[order], female[order]
    depth = np.exp(rng.uniform(np.log(8000), np.log(40000), N)).astype(int)
    sample_ids = [f"S{j + 1:02d}" for j in range(N)]

    # ---- taxon roles -------------------------------------------------------
    share = np.exp(rng.normal(0.0, 1.6, T))
    share /= share.sum()
    idx = rng.permutation(T)
    allzero = idx[:N_ALLZERO]
    nearzero = idx[N_ALLZERO:N_ALLZERO + N_NEARZERO]
    rest = idx[N_ALLZERO + N_NEARZERO:]
    rare = np.r_[allzero, nearzero]

    rank = np.argsort(-share)                          # abundance ranking, high -> low
    eligible = [i for i in rank if i not in set(rare)]
    # bloomers: moderately abundant taxa (ranks 8-60); decliners: abundant (ranks 3-30)
    up_pool = [i for i in eligible[2:35]]
    for _ in range(1000):          # blooms must be big enough to reshape the composition
        up = rng.choice(up_pool, len(UP_FC), replace=False)
        if UP_SHARE[0] <= share[up].sum() <= UP_SHARE[1]:
            break
    else:
        raise RuntimeError("no bloomer set within share window")
    dn_pool = [i for i in eligible[3:30] if i not in set(up)]
    dn = rng.choice(dn_pool, len(DOWN_FC), replace=False)
    truth_idx = np.r_[up, dn]
    nuis_pool = [i for i in eligible[5:120] if i not in set(truth_idx)]
    sex_taxa = rng.choice(nuis_pool, N_SEX, replace=False)

    logfc = np.zeros(T)
    logfc[up] = np.log(UP_FC)
    logfc[dn] = np.log(DOWN_FC)
    sexeff = np.zeros(T)
    sexeff[sex_taxa] = rng.choice([-1, 1], N_SEX) * rng.uniform(0.8, 1.4, N_SEX)

    # near-zero taxa: expected ~0.03-0.2 reads per 20k -> a handful of nonzero entries
    share[nearzero] = rng.uniform(1.5e-6, 1.0e-5, len(nearzero))
    share[allzero] = 0.0

    # structural zeros (biological absence) for ~30% of remaining null taxa
    presence = np.ones(T)
    prot = set(truth_idx) | set(sex_taxa) | set(rare)
    for i in rest:
        if i not in prot and rng.random() < 0.30:
            presence[i] = rng.uniform(0.4, 0.9)
    sd = rng.uniform(0.4, 0.8, T)
    sd[truth_idx] = np.minimum(sd[truth_idx], 0.6)      # keep the true changes detectable

    # ---- absolute abundance & sequencing ----------------------------------
    A = np.zeros((T, N))
    for j in range(N):
        la = (np.log(np.maximum(share, 1e-300)) + logfc * case[j] + sexeff * female[j]
              + rng.normal(0, 1, T) * sd)
        a = np.exp(la) * (rng.random(T) < presence)
        a[allzero] = 0.0
        A[:, j] = a * 1e9                              # arbitrary "cells per g" scale
    P = A / A.sum(axis=0)
    counts = np.zeros((T, N), dtype=np.int64)
    for j in range(N):
        q = P[:, j] * rng.gamma(GAMMA_K, 1.0 / GAMMA_K, T)
        q /= q.sum()
        counts[:, j] = rng.multinomial(depth[j], q)

    # ---- taxon ids (unordered w.r.t. abundance/truth) & taxonomy ----------
    perm = rng.permutation(T)                          # original index -> id number
    ids = [f"ASV{perm[i] + 1:03d}" for i in range(T)]
    lin_choice = rng.integers(0, len(LINEAGES), T)
    unclassified = rng.random(T) < 0.12
    tax_rows = []
    for i in range(T):
        ph, cl, orr, fam, gen = LINEAGES[lin_choice[i]]
        if unclassified[i]:
            gen = "unclassified"
        tax_rows.append((ids[i], "Bacteria", ph, cl, orr, fam, gen))

    return dict(counts=counts, case=case, female=female, depth=depth, sample_ids=sample_ids, ids=ids,
                tax_rows=tax_rows, up=up, dn=dn, sex_taxa=sex_taxa, A=A, P=P, share=share, logfc=logfc,
                allzero=allzero, nearzero=nearzero)


def write():
    w = build()
    D.mkdir(exist_ok=True)
    ids, sids = w["ids"], w["sample_ids"]
    order = sorted(range(T), key=lambda i: ids[i])
    with open(D / "counts.tsv", "w", newline="\n") as f:
        f.write("taxon_id\t" + "\t".join(sids) + "\n")
        for i in order:
            f.write(ids[i] + "\t" + "\t".join(str(int(x)) for x in w["counts"][i]) + "\n")
    with open(D / "metadata.tsv", "w", newline="\n") as f:
        f.write("sample\tgroup\tsex\n")
        for j, s in enumerate(sids):
            f.write(f"{s}\t{'case' if w['case'][j] else 'control'}\t{'F' if w['female'][j] else 'M'}\n")
    with open(D / "taxonomy.tsv", "w", newline="\n") as f:
        f.write("taxon_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\n")
        for i in order:
            f.write("\t".join(w["tax_rows"][i]) + "\n")

    truth_up = sorted(ids[i] for i in w["up"])
    truth_dn = sorted(ids[i] for i in w["dn"])
    others = np.setdiff1d(np.arange(T), np.r_[w["up"], w["dn"], w["allzero"], w["nearzero"]])
    A, P = w["A"], w["P"]
    ctl, cs = w["case"] == 0, w["case"] == 1
    rel_shift = float(np.median(P[others][:, cs].mean(1) / np.maximum(P[others][:, ctl].mean(1), 1e-300)))
    truth = {
        "truth_taxa": sorted(truth_up + truth_dn),
        "up": {ids[i]: float(fc) for i, fc in zip(w["up"], UP_FC)},
        "down": {ids[i]: float(fc) for i, fc in zip(w["dn"], DOWN_FC)},
        "sex_nuisance_taxa_not_in_truth": sorted(ids[i] for i in w["sex_taxa"]),
        "all_zero_taxa": sorted(ids[i] for i in w["allzero"]),
        "near_zero_taxa": sorted(ids[i] for i in w["nearzero"]),
        "note": "All other taxa have identical ABSOLUTE abundance in cases and controls; "
                "their relative abundance falls in cases by the factor below.",
        "median_relative_abundance_ratio_case_over_control_of_unchanged_taxa": rel_shift,
        "total_absolute_load_ratio_case_over_control": float(A[:, cs].sum(0).mean() / A[:, ctl].sum(0).mean()),
    }
    with open(HERE / "truth.json", "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(truth, indent=2))

    task = (
        "You are given a 16S rRNA amplicon sequencing survey of the gut microbiome in 20 cases and 20 controls, "
        "processed into an ASV table: data/counts.tsv (raw read counts; rows are taxa, columns are samples), "
        "data/metadata.tsv (sample, group = case/control, sex) and data/taxonomy.tsv (taxonomic assignment of each "
        "taxon). Sequencing depth varies between samples.\n\n"
        "Which taxa differ in abundance between cases and controls? Adjust for sex if you think it is "
        "appropriate, control the false discovery rate at 5% with the Benjamini-Hochberg procedure, and "
        "report the taxa that pass. Use the taxon ids exactly as they appear in the first column of data/counts.tsv "
        "(e.g. ASV001).\n\n"
        "Return EXACTLY:\n\n<EVAL_ANSWER>\n"
        "{\"differential_taxa\": [\"<taxon_id>\", ...]}\n"
        "</EVAL_ANSWER>"
    )
    graders = [{
        "type": "marker_gene_precision_recall",
        "config": {"answer_field": "differential_taxa", "ground_truth_labels": truth["truth_taxa"],
                   "scoring": {"pass_thresholds": {"precision_at_k": PRECISION_MIN, "recall_at_k": RECALL_MIN}}}}]
    ev = {
        "id": "i8_microbiome_da", "replicates": None, "synthetic": True, "task": task, "graders": graders,
        "data_dir": "data",
        "relevant_skills": ["bio-microbiome-differential-abundance", "bio-microbiome-amplicon-processing",
                            "bio-microbiome-diversity-analysis"],
        "design_source": "Compositionality: a bloom of a few taxa makes all others fall in relative abundance, "
                         "so relative-abundance tests give false positives (Gloor et al. 2017 Front Microbiol; "
                         "Nearing et al. 2022 Nat Commun).",
    }
    with open(HERE / "eval.json", "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(ev, indent=2))
    return truth


if __name__ == "__main__":
    t = write()
    print("truth:", t["truth_taxa"])
    print("median rel-abundance ratio (unchanged taxa, case/control):",
          round(t["median_relative_abundance_ratio_case_over_control_of_unchanged_taxa"], 3))
    print("absolute load ratio:", round(t["total_absolute_load_ratio_case_over_control"], 3))
