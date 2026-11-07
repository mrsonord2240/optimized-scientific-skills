"""Synthetic replica of variantbench `fine_mapping_effector_gene_prioritization`.

Everything here is FABRICATED. Gene symbols are invented (SYN*) so the answer
cannot be recalled from the literature; the original's answer depended partly on
recalling that the effector gene is a negative regulator of interferon signalling,
so we state that in a `gene_description` column instead.

Design source: the benchmark's public `notes` field (what it tests, what the
traps are), not our Skills. Deterministic: seed below.

Run:  python generate.py     -> data/*.tsv and eval.json
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
rng = np.random.default_rng(20260919)

CHR = "chr22"
LOCUS_START, LOCUS_END = 30_000_000, 30_500_000

# --- Genes ---------------------------------------------------------------
DESC = {
    "SYNU1": "negative feedback regulator of type I interferon receptor signalling",
    "SYNU2": "solute carrier family membrane transporter",
    "SYNU3": "DEAD-box RNA helicase",
}

# --- Table 2: 95% credible set --------------------------------------------
# (id, pos, ref, alt, pip, beta, consequence, gene, sift, polyphen)
CS = [
    ("chr22:30231877:A:G", 30231877, "A", "G", 0.36, 0.181, "intron_variant",            "SYNU1", "NA",             "NA"),
    ("chr22:30233402:C:T", 30233402, "C", "T", 0.22, 0.176, "missense_variant",         "SYNU1", "tolerated(0.42)", "benign(0.003)"),
    ("chr22:30229310:G:A", 30229310, "G", "A", 0.12, 0.170, "synonymous_variant",       "SYNU1", "NA",             "NA"),
    ("chr22:30236054:T:C", 30236054, "T", "C", 0.10, 0.165, "intron_variant",           "SYNU1", "NA",             "NA"),
    ("chr22:30201988:A:C", 30201988, "A", "C", 0.08, 0.160, "intergenic_variant",       "-",     "NA",             "NA"),
    ("chr22:30241730:G:T", 30241730, "G", "T", 0.06, 0.158, "3_prime_UTR_variant",      "SYNU1", "NA",             "NA"),
    ("chr22:30188455:C:T", 30188455, "C", "T", 0.05, 0.150, "upstream_gene_variant",    "SYNU2", "NA",             "NA"),
]
PROTEIN_ALTERING = {"missense_variant", "stop_gained", "stop_lost", "start_lost",
                    "frameshift_variant", "inframe_insertion", "inframe_deletion"}
n_pa = sum(1 for v in CS if v[6] in PROTEIN_ALTERING)
assert n_pa == 1
assert sum(v[4] for v in CS) >= 0.95

with open(DATA / "credible_set_variants.tsv", "w", newline="\n") as f:
    f.write("variant_id\tpos\tref\talt\tpip\tbeta\tconsequence\tgene\tgene_description\tsift\tpolyphen\n")
    for vid, pos, ref, alt, pip, beta, cons, gene, sift, pph in CS:
        f.write(f"{vid}\t{pos}\t{ref}\t{alt}\t{pip}\t{beta}\t{cons}\t{gene}\t{DESC.get(gene, '-')}\t{sift}\t{pph}\n")

# --- Table 3: cis-eQTL of SYNU1 expression across the region --------------
# No signal: SYNU1 has no detectable cis-eQTL, credible-set variants included.
cs_pos = {v[1]: v[0] for v in CS}
pos = np.unique(np.concatenate([
    np.sort(rng.integers(LOCUS_START, LOCUS_END, 110)),
    np.array(list(cs_pos)),
]))
with open(DATA / "cis_eqtl_SYNU1.tsv", "w", newline="\n") as f:
    f.write("variant_id\tpos\tbeta\tse\tp\tin_credible_set\n")
    for p in pos:
        se = rng.uniform(0.045, 0.07)
        while True:
            beta = rng.normal(0, se)
            z = beta / se
            from math import erfc, sqrt
            pv = erfc(abs(z) / sqrt(2))
            if pv > 0.02:          # keep the whole region flat: no eQTL of any kind
                break
        vid = cs_pos.get(int(p), f"chr22:{int(p)}")
        f.write(f"{vid}\t{int(p)}\t{beta:.4f}\t{se:.4f}\t{pv:.4g}\t{'yes' if int(p) in cs_pos else 'no'}\n")

# --- Table 1: genes trans-regulated by the locus --------------------------
# 30 interferon-stimulated-like genes, all up with the risk allele; a few unrelated
# genes with weak, non-significant effects. The strongest trans gene is a downstream
# target, not the effector (option F trap).
rows = []
for i in range(1, 31):
    beta = 0.95 - 0.028 * (i - 1) + rng.normal(0, 0.02)
    se = rng.uniform(0.05, 0.09)
    z = beta / se
    from math import erfc, sqrt
    rows.append((f"SYNI{i:02d}", f"chr{rng.integers(1, 22)}", beta, se, erfc(abs(z) / sqrt(2)),
                 "interferon-stimulated gene (type I interferon response)"))
for i in range(1, 11):
    se = rng.uniform(0.06, 0.1)
    beta = rng.normal(0, 0.05)
    z = beta / se
    rows.append((f"SYNX{i:02d}", f"chr{rng.integers(1, 22)}", beta, se, erfc(abs(z) / sqrt(2)),
                 rng.choice(["ribosomal protein", "metabolic enzyme", "cytoskeletal", "unannotated"])))
rows.sort(key=lambda r: -abs(r[2]))
with open(DATA / "trans_genes.tsv", "w", newline="\n") as f:
    f.write("gene\tchr\tbeta\tse\tp\tannotation\n")
    for g, c, b, s, p, a in rows:
        f.write(f"{g}\t{c}\t{b:.4f}\t{s:.4f}\t{p:.3g}\t{a}\n")

# --- eval.json -----------------------------------------------------------
TASK = """You are given three tables for a chr22 disease-associated locus (data/trans_genes.tsv, data/credible_set_variants.tsv, data/cis_eqtl_SYNU1.tsv): (1) the genes trans-regulated by the locus with their effect estimates; (2) the variants in the 95% credible set for the locus signal, each annotated with its functional consequence, the gene it falls in, and (where coding) SIFT and PolyPhen predictions; (3) the cis association of the gene physically spanning the credible set with its own expression across the surrounding region.

Nominate the causal (effector) gene at this locus and the variant class most likely to mediate its effect.

Report the number of protein-altering variants in the credible set, and your interpretation.

Return EXACTLY:

<EVAL_ANSWER>
{"n_protein_altering_variants_in_credible_set": <int>, "answer": "<letter>"}
</EVAL_ANSWER>

where <letter> is one of A, B, C, D, E, F:
A. The effector is the gene spanning the credible set, acting through the intronic lead variant's regulation of that gene's own expression (a cis-expression mechanism).
B. The effector is the gene spanning the credible set, acting through a protein-altering (missense) variant in the credible set, even though that variant is computationally predicted benign and the gene has no detectable cis-eQTL.
C. The effector is whichever gene at the locus has the strongest cis-eQTL, following standard locus-to-gene expression prioritisation.
D. No coding causal mechanism is supported, because protein-altering variants in the credible set are predicted benign, so the effector must act non-coding.
E. The effector is the gene physically closest to the lead variant, by proximity.
F. The effector is the most strongly trans-associated gene, since it carries the strongest association at the locus."""

EVAL = {
    "id": "synthetic_fine_mapping_effector_gene",
    "replicates": "variantbench/fine_mapping_effector_gene_prioritization",
    "synthetic": True,
    "task": TASK,
    "graders": [
        {"type": "numeric_tolerance", "config": {
            "tolerances": {"n_protein_altering_variants_in_credible_set": {"type": "absolute", "value": 0}},
            "ground_truth": {"n_protein_altering_variants_in_credible_set": 1}}},
        {"type": "multiple_choice", "config": {"correct_answer": "B"}},
    ],
    "data_dir": "data",
}
(HERE / "eval.json").write_text(json.dumps(EVAL, indent=2), encoding="utf-8")
print("wrote", sorted(p.name for p in DATA.iterdir()), "+ eval.json")
