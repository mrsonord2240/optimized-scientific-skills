"""Assess how many variants in a VCF require normalization before running bcftools norm.

Usage: python check_normalization.py input.vcf.gz

Counts multiallelic sites and complex (MNP) variants with cyvcf2. Does not detect
indels requiring left-alignment, since that requires reference context -- the count
is a lower bound.
"""
import sys

from cyvcf2 import VCF


def needs_normalization(variant):
    if len(variant.ALT) > 1:
        return True
    ref, alt = variant.REF, variant.ALT[0]
    if len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
        return True
    return False


def main(vcf_path):
    total, needs_norm, multiallelic, mnps = 0, 0, 0, 0
    for variant in VCF(vcf_path):
        total += 1
        if len(variant.ALT) > 1:
            multiallelic += 1
        ref, alt = variant.REF, variant.ALT[0]
        if len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
            mnps += 1
        if needs_normalization(variant):
            needs_norm += 1

    print(f'Total variants: {total}')
    print(f'Needing normalization: {needs_norm} ({needs_norm/total*100:.1f}%)')
    print(f'  Multiallelic sites: {multiallelic}')
    print(f'  MNPs: {mnps}')


if __name__ == '__main__':
    main(sys.argv[1])
