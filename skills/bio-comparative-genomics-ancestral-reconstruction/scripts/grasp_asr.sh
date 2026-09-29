#!/usr/bin/env bash
set -euo pipefail

# GRASP CLI contract (21-Mar-2024). Usage:
#   grasp_asr.sh alignment.fasta species.nwk [output-dir] [threads]
alignment=${1:-alignment.fasta}
tree=${2:-species.nwk}
output_dir=${3:-grasp_run}
threads=${4:-4}

if [[ ! -s "$alignment" ]]; then
  printf 'GRASP alignment is missing or empty: %s\n' "$alignment" >&2
  exit 2
fi
if [[ ! -s "$tree" ]]; then
  printf 'GRASP tree is missing or empty: %s\n' "$tree" >&2
  exit 2
fi
if ! [[ "$threads" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Thread count must be a positive integer: %s\n' "$threads" >&2
  exit 2
fi

mkdir -p "$output_dir"
grasp -a "$alignment" -n "$tree" -o "$output_dir" -j -t "$threads" \
  --save-as FASTA TREE ASR

alignment_name=$(basename "$alignment")
prefix=${alignment_name%.*}
ancestors="$output_dir/${prefix}_ancestors.fa"
ancestor_tree="$output_dir/${prefix}_ancestors.nwk"
asr_json="$output_dir/ASR.json"
for artifact in "$ancestors" "$ancestor_tree" "$asr_json"; do
  if [[ ! -s "$artifact" ]]; then
    printf 'GRASP completed without expected non-empty artifact: %s\n' "$artifact" >&2
    exit 3
  fi
done

# Validate that the advertised ancestor sequences, tree, and JSON are parseable.
python - "$ancestors" "$ancestor_tree" "$asr_json" <<'PY'
import json
import sys
from Bio import AlignIO, Phylo

fasta, tree_path, json_path = sys.argv[1:]
alignment = AlignIO.read(fasta, "fasta")
tree = Phylo.read(tree_path, "newick")
with open(json_path, encoding="utf-8") as handle:
    result = json.load(handle)
if len(alignment) == 0 or alignment.get_alignment_length() == 0:
    raise SystemExit("GRASP ancestor FASTA has no aligned sequences")
if len(tree.get_terminals()) == 0 or not isinstance(result, dict) or not result:
    raise SystemExit("GRASP ancestor tree or ASR.json has no usable records")
print(f"validated {len(alignment)} ancestor sequences x {alignment.get_alignment_length()} columns; "
      f"{len(tree.get_terminals())} tree tips; {len(result)} JSON fields")
PY
