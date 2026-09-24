#!/usr/bin/env python3
'''
Self-test for splice_parsers.py against REAL tool output in test_data/ (GRCh38 panel: ClinVar-pathogenic
canonical splice variants, the GLA c.639+919G>A deep-intronic pseudoexon, ClinVar-benign GLA controls).
Run: python test_splice_parsers.py   (from this directory; pandas + numpy only, no genome needed)
Outputs were produced with SpliceAI 1.3.1 (-A grch38 -D 50 -M 0), Pangolin 1.0.2 (-d 50 -m False, canonical
GENCODE v45 DB), MMSplice 2.4.0 (GENCODE v45 basic GTF); see SKILL.md for the commands.
'''
import math
from pathlib import Path
import tempfile

import numpy as np

from splice_parsers import (build_concordance, classify_delta, display_variant_key, parse_mmsplice_csv,
                            parse_pangolin_vcf, parse_spliceai_vcf, read_input_vcf, unscored_report)
from spliceai_clingen_classify import prepare_spliceai_input

D = Path(__file__).parent / 'test_data'


def test_boundaries():
    # published thresholds are inclusive lower bounds; 0.10 is the BP4 ceiling
    cases = {0.00: 'BP4', 0.09: 'BP4', 0.10: 'BP4', 0.11: 'inconclusive', 0.19: 'inconclusive',
             0.20: 'PP3_supporting', 0.49: 'PP3_supporting', 0.50: 'PP3_supporting_prec0.5',
             0.79: 'PP3_supporting_prec0.5', 0.80: 'PP3_supporting_prec0.8', 1.00: 'PP3_supporting_prec0.8'}
    for score, want in cases.items():
        got = classify_delta([score])[0]
        assert got == want, f'{score}: {got} != {want}'
    # scores as parsed from SpliceAI text ("0.20") must land the same way
    assert classify_delta([float('0.20'), float('0.50'), float('0.80')]).tolist() == \
        ['PP3_supporting', 'PP3_supporting_prec0.5', 'PP3_supporting_prec0.8']
    # a missing score is never benign evidence
    assert classify_delta([np.nan])[0] == 'not_scored'


def test_spliceai_panel():
    s = parse_spliceai_vcf(D / 'spliceai_panel_D50_M0.vcf')
    tp53 = s[s.key == '17:7674292:T>C'].iloc[0]
    assert (tp53.DS_AG, tp53.DS_AL, tp53.delta_max, tp53.DP_AL) == (0.85, 1.00, 1.00, -2)   # hand-read from the VCF line
    gla = s[(s.key == 'X:101399747:C>T') & (s.gene == 'GLA')].iloc[0]
    assert (gla.DS_DG, gla.delta_max) == (0.30, 0.30)                                        # deep-intronic pseudoexon
    assert len(s) == 14 and s.key.nunique() == 9                                             # readthrough/overlapping genes add rows


def test_edge_records_are_not_benign():
    inp = read_input_vcf(D / 'edge_grch38.vcf')
    s = parse_spliceai_vcf(D / 'spliceai_edge_D50_M0.vcf')
    n_masked = s[s.gene == 'TP53'].iloc[0]
    assert math.isnan(n_masked.delta_max)                                    # "." scores are NaN, not 0
    assert classify_delta(s.delta_max.tolist())[0] == 'not_scored'
    # the wrong-REF record has no SpliceAI= tag at all (exit code 0): it must be reported
    assert any('usage_guide_wrong_ref' in m for m in unscored_report(D / 'edge_grch38.vcf', s.key, 'SpliceAI'))
    assert len(inp) == 4                                                     # multi-allelic record split into 2 alleles
    ma = s[s.key.str.startswith('17:7674292:T>')]
    assert sorted(ma.key) == ['17:7674292:T>A', '17:7674292:T>C']
    assert display_variant_key('X:1:' + 'A' * 30 + '>T') == 'X:1:' + 'A' * 24 + '...(30 nt)>T'


def test_symbolic_alleles_are_preserved_as_unscored_before_spliceai():
    with tempfile.TemporaryDirectory() as temp:
        source = Path(temp) / 'input.vcf'
        supported = Path(temp) / 'supported.vcf'
        source.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                          '17\t1\tmixed\tA\tC,<DEL>,*\t.\t.\t.\n')
        skipped = prepare_spliceai_input(source, supported)
        assert skipped == {'17:1:A><DEL>', '17:1:A>*'}
        assert '17\t1\tmixed\tA\tC\t.\t.\t.' in supported.read_text()


def test_pangolin_and_mmsplice():
    p = parse_pangolin_vcf(D / 'pangolin_panel_m0.vcf')
    tp53 = p[p.key == '17:7674292:T>C'].iloc[0]
    assert (round(tp53.gain, 2), tp53.gain_pos, round(tp53.loss, 2), tp53.loss_pos) == (0.72, 47, -0.90, -2)
    assert round(tp53.pangolin_score, 2) == -0.90                            # larger magnitude wins, sign kept
    e = parse_pangolin_vcf(D / 'pangolin_edge_m0.vcf')
    assert e.key.tolist() == ['17:7674292:T>C']                              # ALT[0] only; other alleles get no row
    m = parse_mmsplice_csv(D / 'mmsplice_panel.csv')
    assert m.key.is_unique and len(m) == 8                                   # one row per variant; the deep-intronic one has none
    assert round(m.set_index('key').loc['17:7674292:T>C', 'delta_logit_psi'], 2) == -2.57


def test_concordance_panel():
    T = build_concordance(D / 'panel_grch38.vcf', D / 'spliceai_panel_D50_M0.vcf',
                          D / 'pangolin_panel_m0.vcf', D / 'mmsplice_panel.csv')
    assert len(T) == 9                                                        # outer join: nothing silently dropped
    c = T.set_index('id')
    for name in ('TP53_c.673-2A>G', 'DMD_c.9563+1G>A', 'GLA_c.370-1G>A'):
        assert c.loc[name, 'concordance'] == 'all_predict_disruption', name
    for name in ('GLA_rs2071228_benign', 'GLA_rs151195362_benign'):
        assert c.loc[name, 'concordance'] == 'none_predict_disruption', name
    gla = c.loc['GLA_c.639+919G>A']                                           # MMSplice returns no row for it
    assert (gla.n_scored, gla.n_above, gla.concordance) == (2, 2, 'all_predict_disruption')
    assert c.loc['GLA_c.639+919G>A', 'spliceai_label'] == 'PP3_supporting'    # SpliceAI 0.30, never a stronger tier
    assert not c['concordance'].str.contains('pathogenic').any()              # predictions are not pathogenicity calls
    assert c[c.n_scored < 3].index.tolist() == ['GLA_c.639+919G>A']           # a tool that skipped a variant is visible


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('PASS', name)
