#!/usr/bin/env python3
"""Germline hard filters in Python (cyvcf2): the same thresholds as hard_filter.sh, SNPs and indels separately.

Usage:   python cyvcf2_filter.py in.vcf.gz out.vcf [min_dp=10]
         min_dp=10   site DP floor used by this script (a missing DP passes); min_dp=0 turns it off, which makes
                     the result match hard_filter.sh, whose bcftools expression has no DP term
Needs:   cyvcf2 (Linux/macOS or WSL). Checked on cyvcf2 0.34.0.
Rule:    a missing annotation passes (INFO.get returns None), so hom-alt sites with undefined RankSums survive.
         Records that are neither SNP nor indel (MNP, mixed, symbolic) are counted and kept out of the output with
         a message; normalise first (bcftools norm -m- --atomize) if they matter.
"""
import sys

from cyvcf2 import VCF, Writer

if len(sys.argv) < 3:
    sys.exit(__doc__)
src, dst = sys.argv[1], sys.argv[2]
min_dp = 10.0
for arg in sys.argv[3:]:
    if arg.startswith('min_dp='):
        min_dp = float(arg.split('=', 1)[1])
    else:
        sys.exit(f'unknown argument {arg}')

SNP = dict(fs=60.0, sor=3.0, rprs=-8.0, mq=40.0, mqrs=-12.5)
INDEL = dict(fs=200.0, sor=10.0, rprs=-20.0, mq=None, mqrs=None)


def passes(variant, limits):
    def below(value, floor):
        return value is not None and value < floor

    def above(value, ceiling):
        return value is not None and value > ceiling

    info = variant.INFO
    if (variant.QUAL or 0) < 30:
        return False
    if below(info.get('QD'), 2.0) or above(info.get('FS'), limits['fs']) or above(info.get('SOR'), limits['sor']):
        return False
    if below(info.get('ReadPosRankSum'), limits['rprs']):
        return False
    if limits['mq'] is not None and (below(info.get('MQ'), limits['mq']) or below(info.get('MQRankSum'), limits['mqrs'])):
        return False
    if min_dp and below(info.get('DP'), min_dp):
        return False
    return True


vcf = VCF(src)
writer = Writer(dst, vcf)
counts = dict(snp_in=0, snp_kept=0, indel_in=0, indel_kept=0, other=0)
for variant in vcf:
    if variant.is_snp:
        counts['snp_in'] += 1
        if passes(variant, SNP):
            counts['snp_kept'] += 1
            writer.write_record(variant)
    elif variant.is_indel:
        counts['indel_in'] += 1
        if passes(variant, INDEL):
            counts['indel_kept'] += 1
            writer.write_record(variant)
    else:
        counts['other'] += 1
writer.close()
vcf.close()
print(counts)
if counts['other']:
    print(f"WARNING: {counts['other']} records were neither SNP nor indel and were not written", file=sys.stderr)
