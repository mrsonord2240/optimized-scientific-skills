"""Assess how many variants in a VCF require normalization before running bcftools norm.

Usage: python check_normalization.py input.vcf.gz

Counts multiallelic sites and complex (MNP) variants with cyvcf2. Does not detect
indels requiring left-alignment, since that requires reference context -- the count
is a lower bound.
"""
import sys

from cyvcf2 import VCF


def is_mnp(variant):
    if not variant.ALT:
        return False
    ref, alt = variant.REF, variant.ALT[0]
    return len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt)


def main(vcf_path):
    total, needs_norm, multiallelic, mnps = 0, 0, 0, 0
    for variant in VCF(vcf_path):
        total += 1
        multi = len(variant.ALT) > 1
        mnp = is_mnp(variant)
        multiallelic += multi
        mnps += mnp
        needs_norm += multi or mnp

    print(f'Total variants: {total}')
    if total:
        print(f'Needing normalization: {needs_norm} ({needs_norm/total*100:.1f}%)')
    print(f'  Multiallelic sites: {multiallelic}')
    print(f'  MNPs: {mnps}')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: check_normalization.py input.vcf.gz', file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])
