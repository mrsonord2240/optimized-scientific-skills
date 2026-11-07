#!/usr/bin/env python3
"""CRISPRi / CRISPRa guide-search windows around a TSS (genomic coordinates).

Purpose: return the (start, end) search window for CRISPRi (Dolcetto) or CRISPRa (Calabrese)
         guides, strand-aware. Import the functions or use the CLI.
Inputs:  --tss (int, TSS coordinate), --strand (+ or -), --mode (crispri | crispra)
Usage:   python scripts/tss_windows.py --mode crispri --tss 1000000 --strand +
"""
import argparse

def crispri_window(tss_coord, strand='+'):
    '''Dolcetto convention: search -50 to +300 around the FANTOM5 highest-rank CAGE peak.
    Reason: Sanson 2018 found +25 to +75 nt downstream of the TSS optimal for CRISPRi,
    so rank candidates toward that band; the search is relaxed outward to fill the
    per-gene guide quota when poorly-annotated TSSs leave too few candidates.'''
    if strand == '+':
        return (tss_coord - 50, tss_coord + 300)
    return (tss_coord - 300, tss_coord + 50)

def crispra_window(tss_coord, strand='+'):
    '''Calabrese convention: -150 to -75 upstream of TSS.
    Reason: dCas9-VP64 (and SAM, SunTag) activate maximally when bound
    just upstream of Pol II loading. Horlbeck v2 CRISPRa uses -550 to -25
    (broader, lower per-guide signal). For SAM, prefer Calabrese tightness;
    for SunTag, Horlbeck width is acceptable.
    Caveat: at only 75bp wide, this window routinely fails to contain a full
    6-guide quota's worth of PAM sites passing the GC/poly-T filter -- budget
    for shortfalls (report actual count per gene rather than padding with
    out-of-window guides) or widen to Horlbeck v2 when the quota must be met.'''
    if strand == '+':
        return (tss_coord - 150, tss_coord - 75)
    return (tss_coord + 75, tss_coord + 150)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--mode', choices=['crispri', 'crispra'], required=True)
    ap.add_argument('--tss', type=int, required=True)
    ap.add_argument('--strand', choices=['+', '-'], default='+')
    a = ap.parse_args()
    fn = crispri_window if a.mode == 'crispri' else crispra_window
    lo, hi = fn(a.tss, a.strand)
    print(f'{lo}	{hi}')
