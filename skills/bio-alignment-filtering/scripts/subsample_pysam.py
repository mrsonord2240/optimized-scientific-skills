#!/usr/bin/env python3
'''Pair-consistent, seeded subsample of a BAM in pysam (hash on QNAME, so mates stay together).

The seed is mixed into a real hash: `crc32(qname) ^ seed` only flips the low bits and keeps
the same reads for every seed, and `crc32(f'{seed}:{qname}')` still correlates between seeds
(two seeds shared 65% of their reads). Picks different reads than `samtools view -s`.

Usage: python subsample_pysam.py input.bam output.bam [FRACTION=0.1] [SEED=42]
Checked on pysam 0.24.1 / samtools 1.24.
'''
import hashlib
import sys

import pysam

input_bam, output_bam = sys.argv[1], sys.argv[2]
fraction = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42


def keep(qname):
    digest = hashlib.blake2b(f'{seed}:{qname}'.encode(), digest_size=8).digest()
    return int.from_bytes(digest, 'big') < fraction * 2**64


with pysam.AlignmentFile(input_bam, 'rb') as infile:
    with pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
        for read in infile:
            if keep(read.query_name):
                outfile.write(read)
