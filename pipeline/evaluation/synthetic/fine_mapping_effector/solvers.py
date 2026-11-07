"""Reference solvers for the synthetic fine-mapping task.

`integrated` reads all three tables and should pass. Each `trap_*` applies one
single-evidence heuristic named in the original benchmark's notes and must FAIL.
If any trap passes, the task does not discriminate and the data must be redesigned.
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"
cs = pd.read_csv(D / "credible_set_variants.tsv", sep="\t")
trans = pd.read_csv(D / "trans_genes.tsv", sep="\t")
cis = pd.read_csv(D / "cis_eqtl_SYNU1.tsv", sep="\t")

PA = {"missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
      "inframe_insertion", "inframe_deletion"}
n_pa = int(cs.consequence.isin(PA).sum())
lead = cs.sort_values("pip", ascending=False).iloc[0]
span_gene = lead.gene


def integrated():
    """Coding variant in the spanning gene + flat cis-eQTL + trans program consistent
    with loss of a negative regulator of that program -> B."""
    coding_in_span = cs[cs.consequence.isin(PA) & (cs.gene == span_gene)]
    no_cis = cis[cis.in_credible_set == "yes"].p.min() > 0.05
    span_desc = cs[cs.gene == span_gene].gene_description.iloc[0]
    trans_coherent = (trans[trans.p < 1e-6].annotation.str.contains("interferon").mean() > 0.8
                      and "interferon" in span_desc)
    ans = "B" if (len(coding_in_span) and no_cis and trans_coherent) else "A"
    return {"n_protein_altering_variants_in_credible_set": n_pa, "answer": ans}


def trap_cis_expression():     # A: lead is intronic -> cis-expression mechanism
    return {"n_protein_altering_variants_in_credible_set": n_pa, "answer": "A"}


def trap_strongest_cis_eqtl():  # C
    return {"n_protein_altering_variants_in_credible_set": n_pa, "answer": "C"}


def trap_vep_benign():          # D: predictors say benign -> non-coding
    benign = cs[cs.consequence.isin(PA)].polyphen.str.startswith("benign").all()
    return {"n_protein_altering_variants_in_credible_set": n_pa, "answer": "D" if benign else "B"}


def trap_nearest_gene():        # E
    return {"n_protein_altering_variants_in_credible_set": n_pa, "answer": "E"}


def trap_top_trans():           # F
    return {"n_protein_altering_variants_in_credible_set": n_pa, "answer": "F"}


def trap_count_all_coding():    # right letter, but counts synonymous + UTR as protein-altering
    return {"n_protein_altering_variants_in_credible_set":
            int(cs.consequence.isin(PA | {"synonymous_variant", "3_prime_UTR_variant"}).sum()),
            "answer": "B"}


if __name__ == "__main__":
    ev = HERE / "eval.json"
    solvers = [("integrated (should PASS)", integrated), ("trap A cis-expression", trap_cis_expression),
               ("trap C strongest cis-eQTL", trap_strongest_cis_eqtl), ("trap D VEP-benign", trap_vep_benign),
               ("trap E nearest gene", trap_nearest_gene), ("trap F top trans gene", trap_top_trans),
               ("trap count-all-coding", trap_count_all_coding)]
    ok = True
    for name, fn in solvers:
        ans = fn()
        r = grade(ev, ans)
        per = "".join("P" if g["passed"] else "f" for g in r["graders"])
        print(f"{name:32s} {str(ans):90s} all={r['passed']!s:5s} [{per}]")
        expect = name.startswith("integrated")
        ok &= (r["passed"] == expect)
    print("DISCRIMINATES" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
