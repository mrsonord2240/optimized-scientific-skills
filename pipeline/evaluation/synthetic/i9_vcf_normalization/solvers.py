"""Reference solvers for synthetic task i9_vcf_normalization.

`correct`          : split multi-allelics, trim, left-align against the reference, compare
                     (independent re-implementation, not the code in generate.py).
`bcftools_correct` : same via `bcftools norm -m-any -f` + `bcftools isec` in the bioeval-env
                     docker image (skipped if docker is unavailable).
`trap_*`           : shortcuts that must FAIL.

Usage:  python solvers.py
"""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    out[name] = "".join(buf)
    return out


def read_vcf(path):
    recs = []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        recs.append((f[0], int(f[1]), f[3], f[4].split(",")))
    return recs


REF = read_fasta(D / "reference.fa")
A = read_vcf(D / "caller_A.vcf")
B = read_vcf(D / "caller_B.vcf")


def trim(pos, ref, alt):
    """Remove shared suffix, then shared prefix; keep >=1 base unless an allele becomes empty."""
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def left_align(chrom, pos, ref, alt):
    """Input trimmed (pos, ref, alt) where one allele may be empty. Shift left while the last
    base of the non-empty allele equals the reference base just before, then add anchor base."""
    seq = REF[chrom]
    if ref and alt:                      # SNV / MNV / complex: nothing to shift
        return pos, ref, alt
    while pos > 1 and (ref or alt)[-1] == seq[pos - 2]:
        base = seq[pos - 2]
        ref = base + ref[:-1] if ref else ref
        alt = base + alt[:-1] if alt else alt
        pos -= 1
    anchor = seq[pos - 2]
    return pos - 1, anchor + ref, anchor + alt


def norm_key(chrom, pos, ref, alt):
    p, r, a = trim(pos, ref, alt)
    return (chrom,) + left_align(chrom, p, r, a)


def compare(sa, sb):
    return {"n_both": len(sa & sb), "n_only_A": len(sa - sb), "n_only_B": len(sb - sa)}


def correct():
    sa = {norm_key(c, p, r, a) for c, p, r, alts in A for a in alts}
    sb = {norm_key(c, p, r, a) for c, p, r, alts in B for a in alts}
    return compare(sa, sb)


# ---- traps -----------------------------------------------------------------
def trap_raw():
    """Compare CHROM/POS/REF/ALT exactly as written (records)."""
    sa = {(c, p, r, ",".join(a)) for c, p, r, a in A}
    sb = {(c, p, r, ",".join(a)) for c, p, r, a in B}
    return compare(sa, sb)


def trap_position_only():
    """Match on CHROM:POS."""
    return compare({(c, p) for c, p, _, _ in A}, {(c, p) for c, p, _, _ in B})


def trap_split_only():
    """Split multi-allelics, compare raw alleles (no trimming, no left-alignment)."""
    sa = {(c, p, r, x) for c, p, r, a in A for x in a}
    sb = {(c, p, r, x) for c, p, r, a in B for x in a}
    return compare(sa, sb)


def trap_split_trim_only():
    """Split and trim padding, but do not left-align."""
    sa = {(c,) + trim(p, r, x) for c, p, r, a in A for x in a}
    sb = {(c,) + trim(p, r, x) for c, p, r, a in B for x in a}
    return compare(sa, sb)


def joint_norm(chrom, pos, alleles):
    """Left-align/trim a whole record jointly (what `bcftools norm -f` does without -m)."""
    seq = REF[chrom]
    al = list(alleles)
    changed = True
    while changed:
        changed = False
        if all(al) and len({a[-1] for a in al}) == 1:
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


def trap_norm_no_split():
    """Left-align and trim but leave multi-allelic records intact; compare record keys."""
    def keys(recs):
        out = set()
        for c, p, r, alts in recs:
            q, al = joint_norm(c, p, [r] + alts)
            out.add((c, q, al[0], ",".join(al[1:])))
        return out
    return compare(keys(A), keys(B))


def trap_first_alt_only():
    """Normalise properly but only look at the first ALT of each record."""
    sa = {norm_key(c, p, r, a[0]) for c, p, r, a in A}
    sb = {norm_key(c, p, r, a[0]) for c, p, r, a in B}
    return compare(sa, sb)


# ---- bcftools in docker -----------------------------------------------------
def bcftools_correct():
    script = (
        "set -e; cd /w; mkdir -p /tmp/n; "
        "for X in A B; do "
        "bcftools norm -m-any -f reference.fa caller_$X.vcf 2>/dev/null | bcftools sort -Oz -o /tmp/n/$X.vcf.gz 2>/dev/null; "
        "bcftools index -f /tmp/n/$X.vcf.gz; done; "
        "bcftools isec -c none -p /tmp/n/i /tmp/n/A.vcf.gz /tmp/n/B.vcf.gz; "
        "for f in 0000 0001 0002; do bcftools view -H /tmp/n/i/$f.vcf | wc -l; done"
    )
    r = subprocess.run(["docker", "run", "--rm", "--network", "none", "-v", f"{D}:/w:ro", "bioeval-env",
                        "bash", "-c", script], capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-500:])
    a_only, b_only, both = [int(x) for x in r.stdout.split()[-3:]]
    return {"n_both": both, "n_only_A": a_only, "n_only_B": b_only}


def bcftools_no_split():
    """bcftools norm -f WITHOUT -m, then isec (record level) -- cross-check for trap_norm_no_split."""
    script = (
        "set -e; cd /w; mkdir -p /tmp/n; "
        "for X in A B; do "
        "bcftools norm -f reference.fa caller_$X.vcf 2>/dev/null | bcftools sort -Oz -o /tmp/n/$X.vcf.gz 2>/dev/null; "
        "bcftools index -f /tmp/n/$X.vcf.gz; done; "
        "bcftools isec -c none -p /tmp/n/i /tmp/n/A.vcf.gz /tmp/n/B.vcf.gz; "
        "for f in 0000 0001 0002; do bcftools view -H /tmp/n/i/$f.vcf | wc -l; done"
    )
    r = subprocess.run(["docker", "run", "--rm", "--network", "none", "-v", f"{D}:/w:ro", "bioeval-env",
                        "bash", "-c", script], capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-500:])
    a_only, b_only, both = [int(x) for x in r.stdout.split()[-3:]]
    return {"n_both": both, "n_only_A": a_only, "n_only_B": b_only}


def _run(name, fn, expect, ev, truth):
    try:
        a = fn()
    except Exception as e:  # docker missing etc.
        print(f"{name:26s} SKIPPED ({str(e)[:80]})")
        return True
    r = grade(ev, a)
    per = "".join("P" if g["passed"] else "f" for g in r["graders"])
    dev = {k: a[k] - truth[k] for k in truth}
    print(f"{name:26s} all={r['passed']!s:5s} [{per}] {a}  delta={dev}")
    return r["passed"] == expect


if __name__ == "__main__":
    ev = HERE / "eval.json"
    truth = json.loads((HERE / "truth.json").read_text())
    truth = {k: truth[k] for k in ("n_both", "n_only_A", "n_only_B")}
    print("truth:", truth)
    ok = True
    for name, fn, expect in [("correct (should PASS)", correct, True),
                             ("correct via bcftools", bcftools_correct, True),
                             ("trap raw record compare", trap_raw, False),
                             ("trap CHROM:POS only", trap_position_only, False),
                             ("trap split, no trim/align", trap_split_only, False),
                             ("trap split+trim, no align", trap_split_trim_only, False),
                             ("trap norm, no split", trap_norm_no_split, False),
                             ("trap first ALT only", trap_first_alt_only, False)]:
        ok &= _run(name, fn, expect, ev, truth)
    try:
        print("info: bcftools norm -f (no -m) + isec records:", bcftools_no_split())
    except Exception as e:
        print("info: bcftools no-split skipped:", str(e)[:80])
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
