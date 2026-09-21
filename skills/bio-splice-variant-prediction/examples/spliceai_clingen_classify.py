#!/usr/bin/env python3
'''
SpliceAI VCF annotation with ClinGen SVI 2023 PP3/BP4 evidence labels (research-use decision support,
not a variant classification or diagnosis).

Workflow:
1. Run SpliceAI on a VCF (default window)
2. Parse delta scores; label PP3/BP4 with inclusive published boundaries; missing scores -> not_scored
3. Re-run variants scoring < 0.20 (including 0.00) with a wider window, keep the higher delta
4. Report input variants SpliceAI skipped (it exits 0 on them)

Usage: python spliceai_clingen_classify.py input.vcf genome.fa --build grch38 [--distance 50 --extended 500]
'''
# Reference: spliceai 1.3.1, pandas 2.3 | Verify API if version differs

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

from splice_parsers import PP3_MIN, classify_delta, parse_spliceai_vcf, read_input_vcf, variant_key


def run_spliceai(input_vcf, output_vcf, genome_fa, build, distance=50, mask=0):
    '''Run SpliceAI on a VCF. `build` must match the FASTA (grch37 or grch38).'''
    cmd = ['spliceai', '-I', str(input_vcf), '-O', str(output_vcf), '-R', str(genome_fa),
           '-A', build, '-D', str(distance), '-M', str(mask)]
    subprocess.run(cmd, check=True)


def per_variant(df):
    '''One row per variant: highest delta over the genes SpliceAI annotated (readthrough genes such as
    RPL36A-HNRNPH2 add rows). Choose the MANE gene yourself when the genes disagree.'''
    best = df.sort_values('delta_max', ascending=False, na_position='last').drop_duplicates('key')
    return best.set_index('key')[['gene', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'delta_max']]


def apply_clingen_svi(df):
    '''ClinGen SVI 2023 PP3/BP4 (Walker 2023 AJHG), supporting weight only. 0.5/0.8 are SpliceAI precision
    tiers (Jaganathan 2019), not ACMG strength upgrades; moderate/strong needs functional PS3/BS3.'''
    df = df.copy()
    df['acmg_evidence'] = classify_delta(df['delta_max'])
    return df


def flag_extend_window_candidates(df, threshold=PP3_MIN):
    '''Pseudoexon-creating deep-intronic variants can score 0.00 at the default window, so flag everything
    below the PP3 threshold (scored 0.00 included). Pre-filter the VCF to the intronic/non-coding variants you
    care about if it is large: the wide re-run costs time.'''
    df = df.copy()
    df['extend_window_candidate'] = df['delta_max'] < threshold
    return df


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input_vcf'); ap.add_argument('genome_fa')
    ap.add_argument('--build', required=True, choices=['grch37', 'grch38'])
    ap.add_argument('--distance', type=int, default=50); ap.add_argument('--extended', type=int, default=500)
    ap.add_argument('--out-prefix', default='spliceai_clingen')
    a = ap.parse_args(argv)
    prefix = Path(a.out_prefix)

    narrow = prefix.with_name(prefix.name + f'_D{a.distance}.vcf')
    run_spliceai(a.input_vcf, narrow, a.genome_fa, a.build, distance=a.distance)
    parsed = parse_spliceai_vcf(narrow)
    df = per_variant(parsed)
    inputs = read_input_vcf(a.input_vcf).drop_duplicates('key').set_index('key')
    df = inputs[['id']].join(df)                      # left join: skipped variants stay, as NaN
    df = flag_extend_window_candidates(df)

    cand = df[df['extend_window_candidate']]
    if not cand.empty:
        # subset the ORIGINAL input VCF (keeps its header) to the candidate records
        keys = set(cand.index)
        cand_vcf = prefix.with_name(prefix.name + '_extend_candidates.vcf')
        with open(a.input_vcf) as fin, open(cand_vcf, 'w') as fout:
            for line in fin:
                c = line.split('\t')
                if line.startswith('#') or any(variant_key(c[0], c[1], c[3], alt) in keys for alt in c[4].split(',')):
                    fout.write(line)
        wide = prefix.with_name(prefix.name + f'_D{a.extended}.vcf')
        run_spliceai(cand_vcf, wide, a.genome_fa, a.build, distance=a.extended)
        w = per_variant(parse_spliceai_vcf(wide))['delta_max'].rename('delta_max_wide')
        df = df.join(w)
        df['delta_max'] = df[['delta_max', 'delta_max_wide']].max(axis=1)   # wider window can only add sites
    df = apply_clingen_svi(df)
    df.to_csv(f'{prefix}_classified.tsv', sep='\t')
    print(df['acmg_evidence'].value_counts())
    for key, r in df[df['acmg_evidence'] == 'not_scored'].iterrows():   # no SpliceAI= tag, or "." scores
        print(f'WARNING SpliceAI: no score for {r["id"]} ({key})', file=sys.stderr)
    return df


if __name__ == '__main__':
    main()
