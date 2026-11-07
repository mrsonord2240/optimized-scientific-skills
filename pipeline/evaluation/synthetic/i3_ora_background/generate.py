"""Generate the synthetic ORA-background task (deterministic, seed below).

World: 20,000 annotated genes, 4,000 detected in the assay. Two "tissue program" gene sets
are 70% detected (vs ~19% base rate), so a DE list drawn from detected genes looks enriched
for them against a genome-wide background but not against the detected background.
Two "true" sets carry genuine extra hits in the DE list. All other sets are nulls.
Writes data/, truth.json, eval.json.
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parent
SEED = 20260919
N_GENES, N_DET, N_DE = 20000, 4000, 300
NAME_TRUE, NAME_TISSUE = ["true_A", "true_B"], ["tissue_A", "tissue_B"]


def build(rng):
    genes = np.array([f"GENE{i:05d}" for i in range(1, N_GENES + 1)])
    idx = np.arange(N_GENES)
    perm = rng.permutation(N_GENES)

    # tissue programs: disjoint, 70% detected
    tis = {}
    pool = list(perm)
    for name, size in zip(NAME_TISSUE, (300, 300)):
        tis[name] = np.array([pool.pop() for _ in range(size)])
    det = set()
    for m in tis.values():
        det.update(rng.choice(m, int(round(0.7 * len(m))), replace=False).tolist())
    prog = set(np.concatenate(list(tis.values())).tolist())
    rest = np.array([i for i in idx if i not in prog])
    det.update(rng.choice(rest, N_DET - len(det), replace=False).tolist())
    D = np.array(sorted(det))
    isdet = np.zeros(N_GENES, bool)
    isdet[D] = True
    nonprog = np.array([i for i in idx if i not in prog])

    def draw(size, frac_det):
        """Members drawn from non-program genes with an exact detected fraction."""
        kd = int(round(frac_det * size))
        a = rng.choice(nonprog[isdet[nonprog]], kd, replace=False)
        b = rng.choice(nonprog[~isdet[nonprog]], size - kd, replace=False)
        return np.concatenate([a, b])

    sets = dict(tis)
    sets["true_A"] = draw(90, 0.40)
    sets["true_B"] = draw(120, 0.40)
    specs = []
    for _ in range(6):   # small sets: filtered out under a detected-genes universe
        specs.append((int(rng.integers(15, 41)), 0.15))
    for _ in range(8):   # partially detected
        specs.append((int(rng.integers(80, 201)), 0.38))
    for _ in range(5):   # under-detected
        specs.append((int(rng.integers(150, 301)), 0.08))
    for _ in range(37):  # ordinary sets at base rate
        specs.append((int(rng.integers(40, 301)), 0.19))
    for k, (size, f) in enumerate(specs):
        for _try in range(50):
            m = draw(size, f)
            kd = int(isdet[m].sum())
            if kd <= 7 or kd >= 13:      # keep away from the 10-member filter boundary
                break
        sets[f"null_{k:02d}"] = m
    return genes, D, isdet, sets


def main():
    rng = np.random.default_rng(SEED)
    genes, D, isdet, sets = build(rng)
    names = list(sets)
    order = rng.permutation(len(names))
    ids = {names[j]: f"GS{i + 1:04d}" for i, j in enumerate(order)}
    kinds = {n: ("true" if n.startswith("true") else "tissue" if n.startswith("tissue") else "null") for n in names}

    tr = {n: np.intersect1d(sets[n], D) for n in NAME_TRUE}
    nonT = [n for n in names if kinds[n] == "null"]

    for attempt in range(20000):
        de = []
        for n, kk in zip(NAME_TRUE, (14, 16)):
            de += rng.choice(tr[n], kk, replace=False).tolist()
        de = np.array(sorted(set(de)))
        left = np.setdiff1d(D, de)
        de = np.concatenate([de, rng.choice(left, N_DE - len(de), replace=False)])
        deset = np.zeros(len(genes), bool)
        deset[de] = True
        rows = {}
        for n, m in sets.items():
            Kd = int(isdet[m].sum())
            kd = int(deset[m].sum())
            Kg = len(m)
            rows[n] = dict(Kd=Kd, kd=kd, Kg=Kg,
                           pd=float(hypergeom.sf(kd - 1, N_DET, Kd, N_DE)),
                           pg=float(hypergeom.sf(kd - 1, N_GENES, Kg, N_DE)))
        ok = True
        for n in nonT:
            r = rows[n]
            ok &= r["pd"] >= (0.02 if r["Kd"] >= 13 else 0.05)
            ok &= r["pg"] >= 0.02
        for n in NAME_TISSUE:
            ok &= 0.2 <= rows[n]["pd"] <= 0.95
            ok &= rows[n]["pg"] <= 2e-5
        # raw-p shortcut must catch >=2 nulls (p in [0.02,0.05) on the detected background)
        nraw = sum(1 for n in nonT if rows[n]["Kd"] >= 13 and rows[n]["pd"] < 0.05)
        ok &= nraw >= 2
        if ok:
            break
    else:
        raise SystemExit("no acceptable DE draw")
    print("accepted DE draw at attempt", attempt, "raw-p nulls:", nraw)

    # symbols absent from the annotation in ~a quarter of sets
    gm = {}
    for n, m in sets.items():
        syms = [str(genes[i]) for i in m]
        if rng.random() < 0.25:
            for _ in range(int(rng.integers(1, 4))):
                syms.append(f"LOC{int(rng.integers(100000, 999999))}")
        rng.shuffle(syms)
        gm[n] = syms

    d = HERE / "data"
    d.mkdir(exist_ok=True)

    def w(p, lines):
        (d / p).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    w("all_annotated_genes.txt", [str(g) for g in genes])
    w("detected_genes.txt", [str(genes[i]) for i in D])
    w("de_genes.txt", sorted(str(genes[i]) for i in de))
    gsorted = sorted(names, key=lambda n: ids[n])
    w("pathways.gmt", [f"{ids[n]}\tcurated gene set\t" + "\t".join(gm[n]) for n in gsorted])

    true_ids = sorted(ids[n] for n in NAME_TRUE)
    truth = {"true_sets": true_ids, "tissue_program_sets": sorted(ids[n] for n in NAME_TISSUE),
             "n_annotated": N_GENES, "n_detected": N_DET, "n_de": N_DE, "seed": SEED,
             "design_counts": {ids[n]: rows[n] for n in gsorted}}
    (HERE / "truth.json").write_text(json.dumps(truth, indent=2), encoding="utf-8", newline="\n")

    task = (
        "I ran a differential expression experiment in one tissue and now want to know which curated gene sets are over-represented "
        "among the differentially expressed genes.\n\n"
        "Files:\n"
        "- data/de_genes.txt: the differentially expressed genes (one gene symbol per line, no header).\n"
        "- data/detected_genes.txt: the genes detected in this experiment (one symbol per line).\n"
        "- data/all_annotated_genes.txt: all annotated genes in the genome annotation used (one symbol per line).\n"
        "- data/pathways.gmt: the gene sets, in GMT format (set id, description, then member symbols; tab-separated). "
        "Some members are symbols that are not in the current annotation.\n\n"
        "Test each gene set for over-representation in the DE list with a one-sided hypergeometric test (equivalently, a one-sided Fisher's exact test), "
        "and control the false discovery rate with Benjamini-Hochberg across all gene sets tested. Only test gene sets that have between 10 and 500 members "
        "after restricting the set to the genes in your background; do not test the rest and do not count them in the multiple-testing correction. "
        "Report every gene set with BH-adjusted p < 0.05.\n\n"
        "Return EXACTLY:\n\n<EVAL_ANSWER>\n{\"enriched_gene_sets\": []}\n</EVAL_ANSWER>\n\n"
        "where the list holds the ids of the significant gene sets (first column of the GMT)."
    )
    ev = {"id": "i3_ora_background", "replicates": None, "synthetic": True, "task": task,
          "graders": [{"type": "label_set_jaccard", "config": {
              "scoring": {"pass_threshold": 1}, "answer_field": "enriched_gene_sets",
              "ground_truth_labels": true_ids}}],
          "data_dir": "data",
          "relevant_skills": ["bio-pathway-go-enrichment", "bio-pathway-gsea", "bio-pathway-kegg-pathways",
                              "bio-pathway-reactome", "bio-pathway-wikipathways", "bio-pathway-enrichment-visualization"],
          "design_source": "Over-representation analysis with a whole-genome background instead of the genes detectable in the assay gives spurious enrichment of tissue-expressed gene sets (Timmons et al. 2015 Genome Biol; Wijesooriya et al. 2022 PLoS Comput Biol)."}
    (HERE / "eval.json").write_text(json.dumps(ev, indent=2), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
