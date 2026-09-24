#!/usr/bin/env python3
"""Build the chip-synthesis oligo for one spacer: subpool forward primer + spacer + scaffold.

Purpose: construct the pooled-synthesis oligo for LentiGuide-Puro / LentiCRISPRv2 cloning.
Inputs:  spacer (20-nt string, positional), --subpool (1 or 2; default subpool 1 primer)
Usage:   python scripts/build_oligo.py GGATGGAGACGCATGATTCA --subpool 1
"""
import argparse

def build_oligo(spacer, vector='lentiGuide-Puro', subpool_idx=None):
    '''Construct final oligo for pooled synthesis.

    LentiGuide-Puro / LentiCRISPRv2 use BsmBI (Esp3I) with these overhangs:
        forward: 5'-CACCG[spacer]-3'
        reverse: 5'-AAAC[revcomp(spacer)]C-3'
    For chip synthesis, the spacer is flanked by subpool-specific PCR primers.'''
    subpool_fwd = {
        1: 'GGAAAGGACGAAACACCG',   # subpool 1 forward primer + BsmBI overhang
        2: 'GAGGCACTGGGCAGGTACCG',
    }.get(subpool_idx, 'GGAAAGGACGAAACACCG')
    # First 33 nt of the Chen 2013 sgRNA(F+E) optimized scaffold. NOTE: lentiGuide-Puro (#52963)
    # and lentiCRISPRv2 (#52961) carry the ORIGINAL scaffold; F+E belongs to lentiCRISPRv2-Opti (#163126).
    scaffold_short = 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG'
    oligo = subpool_fwd + spacer + scaffold_short
    if len(oligo) > 200:
        raise ValueError(f'Oligo length {len(oligo)} exceeds the 200 nt design budget; check the vendor limit')
    return oligo

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('spacer')
    ap.add_argument('--subpool', type=int, default=None)
    ap.add_argument('--vector', default='lentiGuide-Puro')
    a = ap.parse_args()
    print(build_oligo(a.spacer.upper(), vector=a.vector, subpool_idx=a.subpool))
