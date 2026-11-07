#!/usr/bin/env python3
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'examples'))
from msa_utils import load_alignment
from henikoff_weights import henikoff_weights
from neff import neff

if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: weights.py alignment output.json')
    alignment = load_alignment(sys.argv[1])
    weights = henikoff_weights(alignment)
    result = {'weights': [{'id': r.id, 'weight': float(w)} for r, w in zip(alignment, weights)],
              'neff_identity_0.62': float(neff(alignment))}
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
