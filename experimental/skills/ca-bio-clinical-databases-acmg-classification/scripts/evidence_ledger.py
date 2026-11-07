"""Non-diagnostic ledger for pre-reviewed, non-patient example criteria."""
import argparse
import json
from pathlib import Path
from acmg_classify import classify_with_subsumption

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("evidence_json")
parser.add_argument("output_json")
args = parser.parse_args()
data = json.loads(Path(args.evidence_json).read_text(encoding="utf-8"))
result = classify_with_subsumption(data["criteria"])
Path(args.output_json).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
