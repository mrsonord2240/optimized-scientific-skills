#!/usr/bin/env python3
'''
Parsers, ClinGen SVI PP3/BP4 labels and cross-tool concordance for SpliceAI, Pangolin and MMSplice output.

Research-use decision support only: labels are computational-evidence tags for expert review, not a
diagnosis or a variant classification. Every un-scored record stays visible as `not_scored`.

Checked on SpliceAI 1.3.1, Pangolin 1.0.2, MMSplice 2.4.0, pandas 2.3.
'''
import re
import numpy as np
import pandas as pd

# ClinGen SVI 2023 (Walker et al.): SpliceAI delta >= 0.20 -> PP3, <= 0.10 -> BP4, both at SUPPORTING weight.
# 0.50 / 0.80 are SpliceAI's own precision tiers (Jaganathan 2019), not ACMG strength upgrades.
BP4_MAX = 0.10
PP3_MIN = 0.20
TIER_50 = 0.50
TIER_80 = 0.80
PANGOLIN_MIN = 0.20      # |score|, Skill convention
MMSPLICE_MIN = 1.0       # |delta_logit_psi|, Skill convention


def norm_chrom(chrom):
    return chrom[3:] if chrom.lower().startswith('chr') else chrom


def variant_key(chrom, pos, ref, alt):
    '''One key for all tools: chromosome without "chr", POS, REF>ALT.'''
    return f'{norm_chrom(str(chrom))}:{int(pos)}:{ref}>{alt}'


def display_variant_key(key, allele_limit=24):
    '''Keep machine keys lossless, but do not put an entire long REF allele in warnings.'''
    chrom, pos, alleles = key.split(':', 2)
    ref, alt = alleles.split('>', 1)
    def shorten(allele):
        return allele if len(allele) <= allele_limit else f'{allele[:allele_limit]}...({len(allele)} nt)'
    return f'{chrom}:{pos}:{shorten(ref)}>{shorten(alt)}'


def classify_delta(delta):
    '''ClinGen SVI PP3/BP4 label per SpliceAI delta_max. Boundaries are inclusive as published
    (0.10 -> BP4, 0.20 -> PP3, 0.50 and 0.80 -> their tiers). NaN (no score) -> not_scored, never BP4.'''
    d = np.round(pd.to_numeric(pd.Series(delta), errors='coerce').to_numpy(dtype=float), 2)
    return np.select(
        [np.isnan(d), d >= TIER_80, d >= TIER_50, d >= PP3_MIN, d <= BP4_MAX],
        ['not_scored', 'PP3_supporting_prec0.8', 'PP3_supporting_prec0.5', 'PP3_supporting', 'BP4'],
        default='inconclusive')


def _num(text):
    return np.nan if text in ('.', '', None) else float(text)


def read_input_vcf(path):
    '''Input variants, one row per ALT allele (multi-allelic records are split).'''
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            for alt in c[4].split(','):
                rows.append({'key': variant_key(c[0], c[1], c[3], alt), 'chrom': c[0], 'pos': int(c[1]),
                             'ref': c[3], 'alt': alt, 'id': c[2]})
    return pd.DataFrame(rows)


def parse_spliceai_vcf(vcf_path):
    '''One row per variant x ALT x gene. "." scores stay NaN. Records without a SpliceAI= tag
    (REF mismatch, deletion longer than 2*D, no gene, unsupported ALT) give no row: find them with
    unscored_report() against the input VCF, SpliceAI exits 0 on them.'''
    rows = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            m = re.search(r'SpliceAI=([^;\t]+)', c[7])
            if not m:
                continue
            for ann in m.group(1).split(','):
                p = ann.split('|')
                if len(p) < 10:
                    continue
                ds = [_num(x) for x in p[2:6]]
                rows.append({'key': variant_key(c[0], c[1], c[3], p[0]), 'gene': p[1],
                             'DS_AG': ds[0], 'DS_AL': ds[1], 'DS_DG': ds[2], 'DS_DL': ds[3],
                             'DP_AG': _num(p[6]), 'DP_AL': _num(p[7]), 'DP_DG': _num(p[8]), 'DP_DL': _num(p[9]),
                             'delta_max': np.nan if all(np.isnan(ds)) else np.nanmax(ds)})
    return pd.DataFrame(rows, columns=['key', 'gene', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL',
                                       'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'delta_max'])


def parse_pangolin_vcf(vcf_path):
    '''One row per variant x gene (Ensembl gene ID). Pangolin scores only ALT[0] of a record and writes
    the max over the four tissue models: gain = largest increase, loss = largest decrease (signed).
    `pangolin_score` is whichever of the two has the larger magnitude.'''
    rows = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            m = re.search(r'Pangolin=([^;\t]+)', c[7])
            if not m:
                continue
            key = variant_key(c[0], c[1], c[3], c[4].split(',')[0])
            for g in m.group(1).split(','):
                p = g.split('|')
                sites = [(int(a), float(b)) for a, b in (x.split(':') for x in p[1:] if re.fullmatch(r'-?\d+:-?[\d.]+(e-?\d+)?', x))]
                gain = max((s for s in sites if s[1] > 0), key=lambda s: s[1], default=(np.nan, 0.0))
                loss = min((s for s in sites if s[1] < 0), key=lambda s: s[1], default=(np.nan, 0.0))
                best = gain if abs(gain[1]) >= abs(loss[1]) else loss
                rows.append({'key': key, 'gene_id': p[0], 'gain_pos': gain[0], 'gain': gain[1],
                             'loss_pos': loss[0], 'loss': loss[1], 'pangolin_score': best[1]})
    return pd.DataFrame(rows, columns=['key', 'gene_id', 'gain_pos', 'gain', 'loss_pos', 'loss', 'pangolin_score'])


def parse_mmsplice_csv(csv_path):
    '''MMSplice writes one row per variant x exon x transcript with ID "chrom:pos:ref>alt".
    Reduce to the row with the largest |delta_logit_psi| per variant.'''
    df = pd.read_csv(csv_path)
    df['key'] = [variant_key(*re.match(r'([^:]+):(\d+):(.+)>(.+)', i).groups()) for i in df['ID']]
    top = df.loc[df['delta_logit_psi'].abs().groupby(df['key']).idxmax()]
    return top[['key', 'gene_name', 'delta_logit_psi', 'pathogenicity']].reset_index(drop=True)


def unscored_report(input_vcf, tool_keys, tool):
    '''Input variants a tool returned nothing for (skipped records exit 0 with no INFO tag).'''
    missing = read_input_vcf(input_vcf).loc[lambda d: ~d['key'].isin(set(tool_keys))]
    return [f'{tool}: no score for {r.id} ({display_variant_key(r.key)})' for r in missing.itertuples()]


def build_concordance(input_vcf, spliceai_vcf=None, pangolin_vcf=None, mmsplice_csv=None):
    '''One row per input variant (outer join: a variant a tool skipped stays, with NaN).
    Above-threshold count is over the tools that scored; fewer than 2 scored -> insufficient_tools.'''
    T = read_input_vcf(input_vcf)[['key', 'id']].drop_duplicates('key').set_index('key')
    if spliceai_vcf:
        s = parse_spliceai_vcf(spliceai_vcf)
        T = T.join(s.groupby('key')['delta_max'].max().rename('spliceai_delta'))   # max over genes; NaN if all "."
    if pangolin_vcf:
        p = parse_pangolin_vcf(pangolin_vcf)
        T = T.join(p.loc[p['pangolin_score'].abs().groupby(p['key']).idxmax()].set_index('key')['pangolin_score'])
    if mmsplice_csv:
        T = T.join(parse_mmsplice_csv(mmsplice_csv).set_index('key')[['delta_logit_psi', 'pathogenicity']])
    thr = {'spliceai_delta': lambda x: x >= PP3_MIN, 'pangolin_score': lambda x: x.abs() >= PANGOLIN_MIN,
           'delta_logit_psi': lambda x: x.abs() >= MMSPLICE_MIN}
    cols = [c for c in thr if c in T]
    T['n_scored'] = sum(T[c].notna().astype(int) for c in cols)
    T['n_above'] = sum((thr[c](T[c]) & T[c].notna()).astype(int) for c in cols)
    T['concordance'] = np.select(
        [T.n_scored < 2, T.n_above == 0, T.n_above == T.n_scored, T.n_above * 2 > T.n_scored],
        ['insufficient_tools', 'none_predict_disruption', 'all_predict_disruption', 'majority_predict_disruption'],
        default='discordant')
    if 'spliceai_delta' in T:
        T['spliceai_label'] = classify_delta(T['spliceai_delta'])
    return T
