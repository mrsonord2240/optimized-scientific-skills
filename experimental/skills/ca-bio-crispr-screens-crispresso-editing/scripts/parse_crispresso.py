#!/usr/bin/env python3
"""Extract key metrics from one CRISPResso2 output directory (CRISPResso_on_<name>/).

Reads CRISPResso_mapping_statistics.txt (7 columns, header + one data row, no percentage column),
CRISPResso_quantification_of_editing_frequency.txt and, when present, CRISPResso2_info.json.
Checked on CRISPResso2 2.3.4, pandas 2.2+.

Usage:
    python parse_crispresso.py <output_dir>     # print the metrics as JSON
    python parse_crispresso.py --selftest       # assert the parser against a synthetic directory
    from parse_crispresso import parse_crispresso   # import from Python (scripts/ on sys.path)

Returned keys: reads_in_input, reads_aligned, mapping_pct, editing_quant ({column: {amplicon: value}}),
and info (parsed CRISPResso2_info.json) when that file exists.
"""
import json
import sys
import tempfile
from pathlib import Path

import pandas as pd


def parse_crispresso(output_dir):
    '''Extract key metrics from CRISPResso output directory.'''
    out = {}
    # Mapping statistics: 7-column, 2-row TSV (header + one data row); no percentage column,
    # so compute mapping_pct from READS ALIGNED / READS IN INPUTS.
    map_stats = pd.read_csv(Path(output_dir) / 'CRISPResso_mapping_statistics.txt', sep='\t').iloc[0]
    out['reads_in_input'] = int(map_stats['READS IN INPUTS'])
    out['reads_aligned'] = int(map_stats['READS ALIGNED'])
    out['mapping_pct'] = out['reads_aligned'] / out['reads_in_input'] * 100
    # Editing quantification
    quant = pd.read_csv(Path(output_dir) / 'CRISPResso_quantification_of_editing_frequency.txt', sep='\t')
    out['editing_quant'] = quant.set_index('Amplicon').to_dict()
    # JSON metadata
    info_path = Path(output_dir) / 'CRISPResso2_info.json'
    if info_path.exists():
        out['info'] = json.loads(info_path.read_text())
    return out


# Real CRISPResso2 2.3.4 output for tests/FANC.Cas9.fastq (250 reads, FANCF amplicon).
_MAP = ('READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\t'
        'N_COMPUTED_NOTALN\tN_CACHED_NOTALN\n250\t250\t235\t198\t37\t15\t0\n')
_QUANT = ('Amplicon\tUnmodified%\tModified%\tReads_in_input\tReads_aligned_all_amplicons\tReads_aligned\t'
          'Unmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\n'
          'Reference\t73.61702128\t26.38297872\t250\t235\t235\t173\t62\t0\t8\t49\t7\n')


def _selftest():
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / 'CRISPResso_mapping_statistics.txt').write_text(_MAP)
        (Path(d) / 'CRISPResso_quantification_of_editing_frequency.txt').write_text(_QUANT)
        r = parse_crispresso(d)
        assert r['reads_in_input'] == 250 and r['reads_aligned'] == 235, r
        assert abs(r['mapping_pct'] - 94.0) < 1e-9, r['mapping_pct']
        assert abs(r['editing_quant']['Modified%']['Reference'] - 26.38297872) < 1e-8, r['editing_quant']
        assert 'info' not in r
        (Path(d) / 'CRISPResso2_info.json').write_text('{"running_info": {"version": "2.3.4"}}')
        assert parse_crispresso(d)['info']['running_info']['version'] == '2.3.4'
    print('selftest OK')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    if sys.argv[1] == '--selftest':
        _selftest()
    else:
        print(json.dumps(parse_crispresso(sys.argv[1]), indent=2, default=str))
