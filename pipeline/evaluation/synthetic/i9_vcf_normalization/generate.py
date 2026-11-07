"""Synthetic eval task i9_vcf_normalization.

Everything here is FABRICATED. Builds a repeat-rich synthetic reference, a pool of
distinct variants with known set membership (both / A-only / B-only), and renders
two VCFs with caller-specific conventions:

  caller A : indels RIGHT-aligned, ~half over-padded with flanking reference bases,
             sites with several alternates written as one multi-allelic record.
  caller B : left-aligned minimal representation, biallelic only.

Truth is by construction (the set each variant was assigned to). Ground truth lives
in truth.json / eval.json only, never in data/.

Run:  python generate.py    -> data/{reference.fa,reference.fa.fai,caller_A.vcf,caller_B.vcf},
                               truth.json, eval.json
"""
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)

SEED = 20260919
rng = random.Random(SEED)

CONTIGS = [("ctg1", 42000), ("ctg2", 31000), ("ctg3", 25000)]
LINE = 60

N_C, N_A, N_B = 204, 52, 41   # designed sizes of both / A-only / B-only (distinct variants)

# ---------------------------------------------------------------- reference ----
def primitive(u):
    n = len(u)
    return not any(n % d == 0 and u == u[:d] * (n // d) for d in range(1, n))


def rand_unit(k):
    while True:
        u = "".join(rng.choice("ACGT") for _ in range(k))
        if primitive(u):
            return u


def rand_seq(n):
    return "".join(rng.choices("ACGT", weights=[3, 2, 2, 3], k=n))


REF = {}      # contig -> str
TRACTS = {}   # contig -> list of (start0, end0, unit)

for name, length in CONTIGS:
    seq, tracts = [], []
    total = 0
    while total < length - 500:
        chunk = rand_seq(rng.randint(90, 330))
        seq.append(chunk)
        total += len(chunk)
        r = rng.random()
        if r < 0.35:      # homopolymer
            b = rng.choices("ATCG", weights=[4, 4, 1, 1])[0]
            unit, n = b, rng.randint(6, 14)
        elif r < 0.65:    # dinucleotide microsatellite
            unit, n = rand_unit(2), rng.randint(5, 12)
        elif r < 0.85:    # trinucleotide
            unit, n = rand_unit(3), rng.randint(4, 8)
        else:             # tetranucleotide
            unit, n = rand_unit(4), rng.randint(3, 7)
        tracts.append((total, total + len(unit) * n, unit))
        seq.append(unit * n)
        total += len(unit) * n
    tail = rand_seq(length - total)
    seq.append(tail)
    s = list(("".join(seq))[:length])
    assert len(s) == length
    # make sure no tract extends into its flanks
    for st, en, u in tracts:
        k = len(u)
        while s[st - 1] == u[-1]:
            s[st - 1] = rng.choice([c for c in "ACGT" if c != u[-1]])
        while s[en] == u[0]:
            s[en] = rng.choice([c for c in "ACGT" if c != u[0]])
    REF[name] = "".join(s)
    TRACTS[name] = tracts

# fasta + fai
fa_lines, fai_lines, offset = [], [], 0
with open(DATA / "reference.fa", "w", newline="\n") as f:
    for name, length in CONTIGS:
        hdr = f">{name}\n"
        f.write(hdr)
        offset += len(hdr)
        fai_lines.append(f"{name}\t{length}\t{offset}\t{LINE}\t{LINE + 1}\n")
        for i in range(0, length, LINE):
            f.write(REF[name][i:i + LINE] + "\n")
        offset += length + (length + LINE - 1) // LINE
with open(DATA / "reference.fa.fai", "w", newline="\n") as f:
    f.writelines(fai_lines)

# ------------------------------------------------------------ normalisation ----
def normalize(contig, pos, alleles):
    """Tan et al. 2015 algorithm. pos is 1-based; alleles = [ref, alt1, ...]."""
    seq = REF[contig]
    al = list(alleles)
    changed = True
    while changed:
        changed = False
        if len({a[-1] for a in al if a}) == 1 and all(al):
            al = [a[:-1] for a in al]
            changed = True
        if any(a == "" for a in al):
            pos -= 1
            al = [seq[pos - 1] + a for a in al]
            changed = True
    while all(len(a) >= 2 for a in al) and len({a[0] for a in al}) == 1:
        al = [a[1:] for a in al]
        pos += 1
    return pos, al


def canon_from_edit(contig, s0, rseg, aseg):
    seq = REF[contig]
    if rseg == "" or aseg == "":
        pos, ref, alt = s0, seq[s0 - 1] + rseg, seq[s0 - 1] + aseg
    else:
        pos, ref, alt = s0 + 1, rseg, aseg
    pos, (ref, alt) = normalize(contig, pos, [ref, alt])
    return pos, ref, alt


def edit_from_canon(pos, ref, alt):
    if len(ref) != len(alt) and ref[0] == alt[0]:
        return pos, ref[1:], alt[1:]
    return pos - 1, ref, alt


def right_align(contig, s, rseg, aseg):
    seq = REF[contig]
    if rseg == "":
        while s < len(seq) and seq[s] == aseg[0]:
            aseg = aseg[1:] + aseg[0]
            s += 1
    elif aseg == "":
        L = len(rseg)
        while s + L < len(seq) and seq[s + L] == seq[s]:
            s += 1
        rseg = seq[s:s + L]
    return s, rseg, aseg


def is_indel(e):
    return e[1] == "" or e[2] == ""


def render_record(contig, edits, padL, padR):
    seq = REF[contig]
    lo = min((s - 1 if is_indel((s, r, a)) else s) for s, r, a in edits) - padL
    hi = max(s + len(r) for s, r, a in edits) + padR
    lo, hi = max(lo, 0), min(hi, len(seq))
    ref = seq[lo:hi]
    alts = [seq[lo:s] + a + seq[s + len(r):hi] for s, r, a in edits]
    assert len(set(alts)) == len(alts) and ref not in alts
    return lo + 1, ref, alts


# --------------------------------------------------------------- placement ----
OCC = {c: [] for c, _ in CONTIGS}
TRACT_ZONES = {c: [(st - 4, en + 4) for st, en, _ in TRACTS[c]] for c, _ in CONTIGS}
USED_TRACTS = set()


def free(contig, lo, hi, zones=False):
    for a, b in OCC[contig]:
        if lo < b and a < hi:
            return False
    if zones:
        for a, b in TRACT_ZONES[contig]:
            if lo < b and a < hi:
                return False
    return True


def pick_contig():
    return rng.choices([c for c, _ in CONTIGS], weights=[l for _, l in CONTIGS])[0]


def other_base(b, avoid=()):
    return rng.choice([c for c in "ACGT" if c != b and c not in avoid])


def pick_tract(min_units=3):
    for _ in range(10000):
        c = pick_contig()
        i = rng.randrange(len(TRACTS[c]))
        st, en, u = TRACTS[c][i]
        if (c, i) in USED_TRACTS or (en - st) // len(u) < min_units:
            continue
        if not free(c, st - 6, en + 6):
            continue
        USED_TRACTS.add((c, i))
        OCC[c].append((st - 6, en + 6))
        return c, st, en, u
    raise RuntimeError("no free tract")


POOL = []   # dicts: contig, edit(s0,rseg,aseg), canon, status, group, padL, padR


def add_variant(contig, edit, status, group=None, pad=None):
    canon = canon_from_edit(contig, *edit)
    v = dict(contig=contig, edit=edit_from_canon(*canon), canon=canon, status=status, group=group, pad=pad)
    POOL.append(v)
    return v


def free_pos(span=1, margin=8, zones=True):
    for _ in range(10000):
        c = pick_contig()
        p = rng.randint(60, len(REF[c]) - 60)
        if free(c, p - margin, p + span + margin, zones=zones):
            OCC[c].append((p - margin, p + span + margin))
            return c, p
    raise RuntimeError("no free position")


def new_snv(status, group=None):
    c, p = free_pos()
    return add_variant(c, (p, REF[c][p], other_base(REF[c][p])), status, group)


def new_mnv(status):
    c, p = free_pos(3)
    n = rng.choice([2, 3])
    ref = REF[c][p:p + n]
    alt = list(ref)
    alt[0] = other_base(ref[0])
    alt[-1] = other_base(ref[-1])
    if n == 3 and rng.random() < 0.5:
        alt[1] = other_base(ref[1])
    return add_variant(c, (p, ref, "".join(alt)), status)


def new_unique_indel(status, kind):
    while True:
        c, p = free_pos(10)
        if kind == "ins":
            e = (p, "", "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6))))
        else:
            L = rng.randint(1, 8)
            e = (p, REF[c][p:p + L], "")
        return add_variant(c, e, status)


def repeat_indel_edit(c, st, en, u, kind):
    k = len(u)
    if k == 1:
        m = rng.choice([1, 1, 2, 3])
    else:
        m = rng.choice([1, 1, 1, 2])
    L = k * m
    if kind == "del":
        L = min(L, ((en - st) // k - 1) * k)
        s = rng.randint(st, en - L)
        return (s, REF[c][s:s + L], "")
    s = rng.randint(st + 1, en - 1)
    return (s, "", REF[c][s:s + k] * m if k > 1 else REF[c][s] * m)


def new_repeat_indel(status, kind, group=None, pad=None):
    c, st, en, u = pick_tract()
    return add_variant(c, repeat_indel_edit(c, st, en, u, kind), status, group, pad)


# ---------------------------------------------------------- multi-allelic ----
# site kinds: snv2 (n alts at one base), insdel (ins+del of one unit at tract start),
# snvindel (SNV on the base before a tract + indel at tract start)
GROUP_SPECS = (
    [("snv2", "CC")] * 4 + [("insdel", "CC")] * 3 + [("snvindel", "CC")] * 3 +
    [("snv2", "CA")] * 2 + [("insdel", "CA")] + [("snvindel", "CA")] +
    [("snv2", "CB")] + [("insdel", "CB")] +
    [("snv2", "AA")] + [("snv3", "CCA")]
)
gid = 0
for kind, sts in GROUP_SPECS:
    gid += 1
    if kind in ("snv2", "snv3"):
        c, p = free_pos()
        ref = REF[c][p]
        alts = []
        for _ in sts:
            alts.append(other_base(ref, avoid=alts))
        for st_, alt in zip(sts, alts):
            add_variant(c, (p, ref, alt), st_, gid)
    else:
        c, st, en, u = pick_tract(min_units=3)
        if kind == "insdel":
            add_variant(c, (st, "", u), sts[0], gid)
            add_variant(c, (st, u, ""), sts[1], gid)
        else:
            base = REF[c][st - 1]
            add_variant(c, (st - 1, base, other_base(base)), sts[0], gid)
            if rng.random() < 0.5:
                add_variant(c, (st, u, ""), sts[1], gid)
            else:
                add_variant(c, (st, "", u), sts[1], gid)

# ---------------------------------------------------------------- decoys ----
# Pairs (one A-only, one B-only) that share the same raw CHROM:POS in the two files
# but are different variants.
for _ in range(6):   # same base, different alt
    v = new_snv("A")
    c, p, ref, alt = v["contig"], v["edit"][0], v["edit"][1], v["edit"][2]
    add_variant(c, (p, ref, other_base(ref, avoid=[alt])), "B")
for _ in range(5):   # A-only right-aligned repeat indel (unpadded); B-only SNV at its raw POS
    kind = rng.choice(["ins", "del"])
    v = new_repeat_indel("A", kind, pad=(0, 0))
    s, r, a = right_align(v["contig"], *v["edit"])
    p0 = s - 1   # raw POS (1-based) == s, i.e. 0-based index s-1 is the anchor base
    c = v["contig"]
    add_variant(c, (p0, REF[c][p0], other_base(REF[c][p0])), "B")
for _ in range(4):   # B-only left-aligned repeat indel; A-only SNV at its raw POS
    kind = rng.choice(["ins", "del"])
    v = new_repeat_indel("B", kind)
    pos = v["canon"][0]
    c = v["contig"]
    add_variant(c, (pos - 1, REF[c][pos - 1], other_base(REF[c][pos - 1])), "A")

# ------------------------------------------------------- remaining singles ----
def count(status):
    return sum(1 for v in POOL if v["status"] == status)


def draw_kind():
    return rng.choices(["snv", "ins", "del"], weights=[60, 16, 24])[0]


def new_single(status, kind):
    if kind == "snv":
        return new_snv(status)
    if kind == "mnv":
        return new_mnv(status)
    if rng.random() < 0.75:
        return new_repeat_indel(status, kind)
    return new_unique_indel(status, kind)


need = {"C": N_C - count("C"), "A": N_A - count("A"), "B": N_B - count("B")}
assert all(x > 0 for x in need.values()), need
n_mnv = {"C": 5, "A": 1, "B": 1}
for status in ("A", "B", "C"):
    kinds = ["mnv"] * n_mnv[status] + [draw_kind() for _ in range(need[status] - n_mnv[status])]
    for k in kinds:
        new_single(status, k)

# ------------------------------------------------------------- uniqueness ----
keys = [(v["contig"],) + v["canon"] for v in POOL]
assert len(set(keys)) == len(keys), "duplicate canonical variants"
for st in "ABC":
    assert sum(1 for v in POOL if v["status"] == st) == {"A": N_A, "B": N_B, "C": N_C}[st], st

# ----------------------------------------------------------------- render ----
CORD = {c: i for i, (c, _) in enumerate(CONTIGS)}


def pad_for(edits, prob_indel=0.45):
    if any(is_indel(e) for e in edits) and len(edits) == 1:
        p = prob_indel
    elif len(edits) > 1:
        p = 0.5
    else:
        p = 0.5 if len(edits[0][1]) > 1 else 0.0   # MNV: sometimes over-padded; SNV never
    if rng.random() >= p:
        return 0, 0
    while True:
        l, r = rng.randint(0, 5), rng.randint(0, 4)
        if l + r > 0:
            return l, r


def genotype(multi=False):
    return "1/2" if multi else rng.choice(["0/1", "0/1", "0/1", "1/1"])


def vcf_header(caller):
    h = ["##fileformat=VCFv4.2", f"##source={caller}"]
    h += [f"##contig=<ID={c},length={l}>" for c, l in CONTIGS]
    h += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1"]
    return h


def emit(records, path, caller):
    records.sort(key=lambda r: (CORD[r[0]], r[1], len(r[2]), r[3]))
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(vcf_header(caller)) + "\n")
        for c, pos, ref, alts, gt in records:
            qual = f"{rng.uniform(20, 900):.1f}"
            dp = rng.randint(8, 90)
            f.write(f"{c}\t{pos}\t.\t{ref}\t{','.join(alts)}\t{qual}\tPASS\t.\tGT:DP\t{gt}:{dp}\n")


# caller A: group multi-allelic sites; right-align indels; over-pad some
recA = []
groups = {}
for v in POOL:
    if "A" in v["status"] or "C" in v["status"]:
        if v["group"] is not None:
            groups.setdefault(v["group"], []).append(v)
        else:
            groups.setdefault(("s", len(groups)), []).append(v)
for key, vs in groups.items():
    c = vs[0]["contig"]
    edits = [right_align(c, *v["edit"]) for v in vs]
    padded = vs[0]["pad"] if len(vs) == 1 and vs[0]["pad"] is not None else pad_for(edits)
    pos, ref, alts = render_record(c, edits, *padded)
    recA.append((c, pos, ref, alts, genotype(len(alts) > 1)))
# caller B: canonical, biallelic
recB = []
for v in POOL:
    if v["status"] in ("B", "C"):
        pos, ref, alt = v["canon"]
        recB.append((v["contig"], pos, ref, [alt], genotype()))

emit(recA, DATA / "caller_A.vcf", "callerA")
emit(recB, DATA / "caller_B.vcf", "callerB")

# -------------------------------------------------- self-check + write truth ----
def norm_set(records):
    out = []
    for c, pos, ref, alts, _ in records:
        for a in alts:
            p, (r, x) = normalize(c, pos, [ref, a])
            out.append((c, p, r, x))
    assert len(set(out)) == len(out)
    return set(out)


SA, SB = norm_set(recA), norm_set(recB)
truth = {"n_both": len(SA & SB), "n_only_A": len(SA - SB), "n_only_B": len(SB - SA)}
assert truth == {"n_both": N_C, "n_only_A": N_A, "n_only_B": N_B}, truth

n_multi = sum(1 for r in recA if len(r[3]) > 1)
raw_pos_shared = len({(r[0], r[1]) for r in recA} & {(r[0], r[1]) for r in recB})
(HERE / "truth.json").write_text(json.dumps({
    "seed": SEED, **truth,
    "n_records_A": len(recA), "n_records_B": len(recB), "n_multiallelic_records_A": n_multi,
    "n_variants_A": len(SA), "n_variants_B": len(SB),
    "shared_raw_chrom_pos": raw_pos_shared,
}, indent=2), encoding="utf-8")

TASK = """Two variant callers were run on the same sample against the same reference, and I have their outputs: data/caller_A.vcf and data/caller_B.vcf (plain-text VCFs, one sample each), with the reference they were called against in data/reference.fa (indexed, data/reference.fa.fai). Every record in both files is a call we want to keep; nothing needs filtering.

I want to know how well the two callers agree. Please count:

- n_both: variants called by both callers
- n_only_A: variants called only by caller A
- n_only_B: variants called only by caller B

Counting conventions:
- Two calls are the same variant if they describe the same change to the reference sequence at the same place, regardless of how each caller happened to write it in the VCF.
- A record with several ALT alleles counts as one variant per ALT allele.
- Genotype, genotype quality, depth and QUAL are ignored; only whether the variant was called matters.
- Multi-base substitutions (MNVs) are counted as single variants; do not break them into separate SNVs.
- Count distinct variants, so n_both + n_only_A is the number of distinct variants in caller A and n_both + n_only_B is the number in caller B.

Return EXACTLY:

<EVAL_ANSWER>
{"n_both": 0, "n_only_A": 0, "n_only_B": 0}
</EVAL_ANSWER>
"""

eval_spec = {
    "id": "i9_vcf_normalization",
    "replicates": None,
    "synthetic": True,
    "task": TASK,
    "graders": [{"type": "numeric_tolerance", "config": {
        "tolerances": {k: {"type": "absolute", "value": 0} for k in truth},
        "ground_truth": truth}}],
    "data_dir": "data",
    "relevant_skills": ["bio-variant-normalization", "bio-vcf-manipulation", "bio-vcf-basics",
                        "bio-vcf-statistics"],
    "design_source": "Callsets from different callers write the same variant differently "
                     "(non-left-aligned/over-padded indels, multi-allelic records); comparison needs "
                     "split + left-align + trim first (Tan et al. 2015 Bioinformatics 31:2202, vt normalize; bcftools norm).",
}
(HERE / "eval.json").write_text(json.dumps(eval_spec, indent=2) + "\n", encoding="utf-8")
print("truth:", truth, "| records A/B:", len(recA), len(recB), "| multi-allelic A records:", n_multi)
