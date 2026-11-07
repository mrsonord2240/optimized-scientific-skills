#!/usr/bin/env python3
"""Per-record VCF statistics with cyvcf2.

Usage: vcf_stats.py <input.vcf.gz> [--per-sample] [--afs]
Always prints counts, PASS/filtered, Ti/Tv (allele-level) and mean QUAL.
--per-sample  also prints het/hom-alt/missing per sample (cyvcf2 gt_types: 0 hom-ref, 1 het, 2 unknown, 3 hom-alt)
--afs         also prints an INFO/AF spectrum (first ALT's AF)
Checked on cyvcf2 0.30+.
"""
import sys

from cyvcf2 import VCF

TI_PAIRS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def main():
    args = sys.argv[1:]
    if not args or args[0].startswith("-"):
        print(__doc__, file=sys.stderr)
        return 1
    path = args[0]
    per_sample = "--per-sample" in args
    afs = "--afs" in args

    vcf = VCF(path)
    samples = vcf.samples
    n = len(samples)
    hom_ref, het, hom_alt, missing = [0] * n, [0] * n, [0] * n, [0] * n
    bins = {"rare(<1%)": 0, "low(1-5%)": 0, "common(5-50%)": 0, "frequent(>50%)": 0}
    s = dict(total=0, snps=0, indels=0, other=0, passed=0, filtered=0,
             ti=0, tv=0, qual_sum=0.0, qual_n=0)

    for v in vcf:
        s["total"] += 1
        if v.is_snp:
            s["snps"] += 1
            # Ti/Tv is allele-level: a multiallelic SNP record counts once per ALT allele.
            for alt in v.ALT:
                if (v.REF, alt) in TI_PAIRS:
                    s["ti"] += 1
                else:
                    s["tv"] += 1
        elif v.is_indel:
            s["indels"] += 1
        else:
            s["other"] += 1
        if v.FILTER is None:
            s["passed"] += 1
        else:
            s["filtered"] += 1
        if v.QUAL is not None:  # QUAL=0 is observed quality, not missing
            s["qual_sum"] += v.QUAL
            s["qual_n"] += 1
        if per_sample:
            for i, gt in enumerate(v.gt_types):
                if gt == 0:
                    hom_ref[i] += 1
                elif gt == 1:
                    het[i] += 1
                elif gt == 3:
                    hom_alt[i] += 1
                else:
                    missing[i] += 1
        if afs:
            af = v.INFO.get("AF")
            if af is not None:
                af = af[0] if isinstance(af, (tuple, list)) else af
                if af < 0.01:
                    bins["rare(<1%)"] += 1
                elif af < 0.05:
                    bins["low(1-5%)"] += 1
                elif af < 0.5:
                    bins["common(5-50%)"] += 1
                else:
                    bins["frequent(>50%)"] += 1
    vcf.close()

    print("=== VCF Statistics ===")
    print(f"Total variants: {s['total']}  SNPs: {s['snps']}  Indels: {s['indels']}  Other: {s['other']}")
    print(f"FILTER null (PASS or '.'): {s['passed']}  Filtered: {s['filtered']}")
    if s["tv"]:
        print(f"Transitions: {s['ti']}  Transversions: {s['tv']}  Ti/Tv: {s['ti'] / s['tv']:.2f}")
    if s["qual_n"]:
        print(f"Mean QUAL: {s['qual_sum'] / s['qual_n']:.1f}")
    if per_sample:
        print("\nsample\thet/hom_alt\thet\thom_alt\tmissing")
        for i, name in enumerate(samples):
            ratio = het[i] / hom_alt[i] if hom_alt[i] else 0
            print(f"{name}\t{ratio:.2f}\t{het[i]}\t{hom_alt[i]}\t{missing[i]}")
    if afs:
        print("\nAF spectrum:")
        for k, c in bins.items():
            print(f"{k}\t{c}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
