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

from splice_parsers import (PP3_MIN, classify_delta, display_variant_key, parse_spliceai_vcf,
                            read_input_vcf, variant_key)


def run_spliceai(input_vcf, output_vcf, genome_fa, build, distance=50, mask=0):
    '''Run SpliceAI on a VCF. `build` must match the FASTA (grch37 or grch38).'''
    cmd = ['spliceai', '-I', str(input_vcf), '-O', str(output_vcf), '-R', str(genome_fa),
           '-A', build, '-D', str(distance), '-M', str(mask)]
    subprocess.run(cmd, check=True)


def prepare_spliceai_input(input_vcf, output_vcf):
    '''Copy a VCF without unsupported ``*`` or symbolic ALT alleles.

    SpliceAI aborts a whole batch on a symbolic ALT.  Retain supported alleles from a
    mixed record, and return the skipped allele keys so they remain visible as
    ``not_scored`` in the final table rather than discarding the batch.
    '''
    skipped = set()
    with open(input_vcf) as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            supported = []
            for alt in fields[4].split(','):
                if alt == '*' or (alt.startswith('<') and alt.endswith('>')):
                    skipped.add(variant_key(fields[0], fields[1], fields[3], alt))
                else:
                    supported.append(alt)
            if supported:
                fields[4] = ','.join(supported)
                fout.write('\t'.join(fields) + '\n')
    return skipped


def per_variant(df):
    '''One row per variant, retaining every SpliceAI gene annotation.

    ``top_score_gene`` is only the annotation with the highest delta; it is not a
    transcript selection. Resolve the disease-relevant/MANE transcript separately.
    '''
    best = df.sort_values('delta_max', ascending=False, na_position='last').drop_duplicates('key')
    best = best.rename(columns={'gene': 'top_score_gene'})
    genes = (df.dropna(subset=['gene']).groupby('key')['gene'].agg(
        lambda names: ','.join(dict.fromkeys(names))).rename('annotated_genes'))
    return best.set_index('key')[['top_score_gene', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'delta_max']].join(genes)


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

    supported = prefix.with_name(prefix.name + '_supported.vcf')
    skipped = prepare_spliceai_input(a.input_vcf, supported)
    narrow = prefix.with_name(prefix.name + f'_D{a.distance}.vcf')
    run_spliceai(supported, narrow, a.genome_fa, a.build, distance=a.distance)
    parsed = parse_spliceai_vcf(narrow)
    df = per_variant(parsed)
    inputs = read_input_vcf(a.input_vcf).drop_duplicates('key').set_index('key')
    df = inputs[['id']].join(df)                      # left join: skipped variants stay, as NaN
    df = flag_extend_window_candidates(df)

    cand = df[df['extend_window_candidate']]
    if not cand.empty:
        # subset the supported input VCF (keeps its header and cannot reintroduce symbolic ALTs)
        keys = set(cand.index)
        cand_vcf = prefix.with_name(prefix.name + '_extend_candidates.vcf')
        with open(supported) as fin, open(cand_vcf, 'w') as fout:
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
        reason = ' (unsupported ALT pre-filtered)' if key in skipped else ''
        print(f'WARNING SpliceAI: no score for {r["id"]} ({display_variant_key(key)}){reason}', file=sys.stderr)
    return df


if __name__ == '__main__':
    main()
