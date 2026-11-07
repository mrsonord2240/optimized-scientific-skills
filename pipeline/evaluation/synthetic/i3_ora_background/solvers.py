"""Reference solvers for the synthetic ORA-background task.

`correct`     : background = genes detected in the experiment.
`trap_*`      : shortcuts that must FAIL (whole-genome background, union-of-pathway-genes
                background, no multiple-testing correction).
`variant_*`   : legitimate implementation variants that must PASS (robustness of the key).

Usage: python solvers.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, hypergeom

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"


def _read(p):
    return [x.strip() for x in (D / p).read_text(encoding="utf-8").splitlines() if x.strip()]


DE = set(_read("de_genes.txt"))
DETECTED = set(_read("detected_genes.txt"))
ANNOTATED = set(_read("all_annotated_genes.txt"))
GMT = {}
for line in (D / "pathways.gmt").read_text(encoding="utf-8").splitlines():
    f = line.split("\t")
    GMT[f[0]] = set(f[2:])
UNION = set().union(*GMT.values())


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def ora(background, lo=10, hi=500, restrict_sets=True, test="hyper"):
    """Return {set_id: (K, k, p, q)} for sets passing the size filter."""
    bg = set(background)
    de = DE & bg
    N, n = len(bg), len(de)
    rows = {}
    for sid, m in GMT.items():
        mm = (m & bg) if restrict_sets else m
        if not (lo <= len(mm) <= hi):
            continue
        K = len(mm)
        k = len(de & mm)
        if test == "hyper":
            p = hypergeom.sf(k - 1, N, K, n)
        else:
            p = fisher_exact([[k, n - k], [K - k, N - n - (K - k)]], alternative="greater")[1]
        rows[sid] = [K, k, p]
    ids = list(rows)
    q = bh([rows[i][2] for i in ids])
    for i, qq in zip(ids, q):
        rows[i].append(float(qq))
    return rows


def sig(rows, alpha=0.05, adjust=True):
    return sorted(s for s, (K, k, p, q) in rows.items() if (q if adjust else p) < alpha)


def correct():
    return {"enriched_gene_sets": sig(ora(DETECTED))}


def trap_whole_genome():
    return {"enriched_gene_sets": sig(ora(ANNOTATED))}


def trap_pathway_union():
    return {"enriched_gene_sets": sig(ora(UNION))}


def trap_no_correction():
    return {"enriched_gene_sets": sig(ora(DETECTED), adjust=False)}


def variant_fisher():
    return {"enriched_gene_sets": sig(ora(DETECTED, test="fisher"))}


def variant_no_size_filter():
    return {"enriched_gene_sets": sig(ora(DETECTED, lo=1, hi=10**9))}


def variant_filter_gt10():
    return {"enriched_gene_sets": sig(ora(DETECTED, lo=11))}


def trap_sets_not_restricted():
    """Set size K taken from the raw GMT list, not restricted to the background (a mistake, must fail)."""
    return {"enriched_gene_sets": sig(ora(DETECTED, restrict_sets=False))}


def variant_filter_on_raw_size():
    """Size filter applied to the raw GMT size, K after restriction."""
    rows = ora(DETECTED, lo=1, hi=10**9)
    keep = {s: r for s, r in rows.items() if 10 <= len(GMT[s]) <= 500}
    ids = list(keep)
    q = bh([keep[i][2] for i in ids])
    return {"enriched_gene_sets": sorted(i for i, qq in zip(ids, q) if qq < 0.05)}


def variant_hits_only():
    """Common shortcut that is still legitimate: test only sets with >=1 DE hit."""
    rows = {s: r for s, r in ora(DETECTED).items() if r[1] >= 1}
    ids = list(rows)
    q = bh([rows[i][2] for i in ids])
    return {"enriched_gene_sets": sorted(i for i, qq in zip(ids, q) if qq < 0.05)}


def variant_bonferroni():
    rows = ora(DETECTED)
    m = len(rows)
    return {"enriched_gene_sets": sorted(s for s, r in rows.items() if r[2] * m < 0.05)}


if __name__ == "__main__":
    import json
    ev = HERE / "eval.json"
    truth = json.loads((HERE / "truth.json").read_text())

    for name, bg in [("detected", DETECTED), ("genome", ANNOTATED), ("union", UNION)]:
        rows = ora(bg)
        s = sig(rows)
        print(f"[{name}] N={len(bg & set(bg))} n_DE_in_bg={len(DE & bg)} tested={len(rows)} significant={len(s)}")
    rows = ora(DETECTED)
    qs = sorted(r[3] for s, r in rows.items())
    print("correct-background q, smallest 6:", ["%.2e" % x for x in qs[:6]])
    print("closest non-significant: smallest q among non-true =",
          "%.3f" % min(r[3] for s, r in rows.items() if s not in truth["true_sets"]),
          "; smallest p among non-true = %.4f" % min(r[2] for s, r in rows.items() if s not in truth["true_sets"]))
    print("largest true-set q = %.2e" % max(rows[s][3] for s in truth["true_sets"]))
    for nm, sid in [("tissue", truth["tissue_program_sets"])]:
        for s in sid:
            print(f"  tissue set {s}: detected bg p={rows[s][2]:.3f}; genome bg p={ora(ANNOTATED)[s][2]:.2e} q={ora(ANNOTATED)[s][3]:.2e}")
    print()

    ok = True
    cases = [("correct (should PASS)", correct, True),
             ("trap whole-genome bg", trap_whole_genome, False),
             ("trap pathway-union bg", trap_pathway_union, False),
             ("trap no FDR correction", trap_no_correction, False),
             ("variant Fisher exact", variant_fisher, True),
             ("variant no size filter", variant_no_size_filter, True),
             ("variant filter >10", variant_filter_gt10, True),
             ("trap sets not restricted", trap_sets_not_restricted, False),
             ("variant filter on raw size", variant_filter_on_raw_size, True),
             ("variant hits-only sets", variant_hits_only, True),
             ("variant Bonferroni", variant_bonferroni, True)]
    for name, fn, expect in cases:
        a = fn()
        r = grade(ev, a)
        per = "".join("P" if g["passed"] else "f" for g in r["graders"])
        print(f"{name:28s} all={r['passed']!s:5s} [{per}] n={len(a['enriched_gene_sets'])} {a['enriched_gene_sets']}")
        ok &= (r["passed"] == expect)
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
