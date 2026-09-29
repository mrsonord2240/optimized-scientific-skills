#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat >&2 <<'EOF'
Usage:
  run_chimeric_eclip.sh --reads FASTQ --hyb-db NAME --run-id ID \
    --output-dir DIR [--expressed-mirnas TSV --expression-threshold NUMBER] \
    [--hyb-bin PATH] [--replicates N] [--replace]

The database must already be installed beneath $HYB_HOME/data/db and must have
a Bowtie2 index. The expression TSV, when supplied, must have this header:
mirna_id<TAB>expression_value<TAB>expression_unit<TAB>expression_source

Outputs are published as one directory only after every Hyb replicate and the
consensus parser succeed. Existing output directories are refused unless
--replace is explicit.
EOF
  exit 64
}

reads=; hyb_db=; run_id=; output_dir=; expressed_mirnas=; expression_threshold=
hyb_bin=hyb; replicates=2; replace=0
while (($#)); do
  case "$1" in
    --reads) [[ $# -ge 2 ]] || usage; reads=$2; shift 2 ;;
    --hyb-db) [[ $# -ge 2 ]] || usage; hyb_db=$2; shift 2 ;;
    --run-id) [[ $# -ge 2 ]] || usage; run_id=$2; shift 2 ;;
    --output-dir) [[ $# -ge 2 ]] || usage; output_dir=$2; shift 2 ;;
    --expressed-mirnas) [[ $# -ge 2 ]] || usage; expressed_mirnas=$2; shift 2 ;;
    --expression-threshold) [[ $# -ge 2 ]] || usage; expression_threshold=$2; shift 2 ;;
    --hyb-bin) [[ $# -ge 2 ]] || usage; hyb_bin=$2; shift 2 ;;
    --replicates) [[ $# -ge 2 ]] || usage; replicates=$2; shift 2 ;;
    --replace) replace=1; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$reads" && -n "$hyb_db" && -n "$run_id" && -n "$output_dir" ]] || usage
[[ -r "$reads" && -s "$reads" ]] || { echo "Unreadable or empty FASTQ: $reads" >&2; exit 66; }
[[ "$hyb_db" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Invalid Hyb database name" >&2; exit 64; }
[[ "$run_id" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Invalid run id" >&2; exit 64; }
[[ "$replicates" =~ ^[0-9]+$ && "$replicates" -ge 2 ]] || {
  echo "--replicates must be an integer >= 2 for the consensus policy" >&2; exit 64;
}
[[ -n "${HYB_HOME:-}" && -d "$HYB_HOME/data/db" ]] || {
  echo "HYB_HOME must identify a current Hyb installation with data/db" >&2; exit 69;
}
[[ -e "$HYB_HOME/data/db/${hyb_db}.1.bt2" || -e "$HYB_HOME/data/db/${hyb_db}.1.bt2l" ]] || {
  echo "No Bowtie2 index found for named Hyb database '$hyb_db'" >&2; exit 66;
}
command -v "$hyb_bin" >/dev/null 2>&1 || { echo "Hyb executable not found: $hyb_bin" >&2; exit 69; }
command -v python3 >/dev/null 2>&1 || { echo "python3 is required" >&2; exit 69; }

if [[ -n "$expressed_mirnas" ]]; then
  [[ -r "$expressed_mirnas" && -s "$expressed_mirnas" ]] || {
    echo "Unreadable or empty expression TSV: $expressed_mirnas" >&2; exit 66;
  }
  [[ -n "$expression_threshold" ]] || {
    echo "--expression-threshold is required with --expressed-mirnas" >&2; exit 64;
  }
elif [[ -n "$expression_threshold" ]]; then
  echo "--expression-threshold requires --expressed-mirnas" >&2; exit 64
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
reads=$(readlink -f -- "$reads")
if [[ -n "$expressed_mirnas" ]]; then expressed_mirnas=$(readlink -f -- "$expressed_mirnas"); fi
output_parent=$(dirname -- "$output_dir"); output_name=$(basename -- "$output_dir")
[[ -n "$output_name" && "$output_name" != "." && "$output_name" != ".." && "$output_name" != "/" ]] || {
  echo "Unsafe output directory name" >&2; exit 64
}
mkdir -p -- "$output_parent"; output_parent=$(cd -- "$output_parent" && pwd -P)
output_dir="$output_parent/$output_name"
if [[ -e "$output_dir" && "$replace" -ne 1 ]]; then
  echo "Output exists; use --replace to replace it atomically: $output_dir" >&2; exit 73
fi

stage=$(mktemp -d -- "$output_parent/.${output_name}.stage.XXXXXXXX"); backup=
cleanup() {
  status=$?
  if [[ -n "$backup" && -e "$backup" && ! -e "$output_dir" ]]; then mv -- "$backup" "$output_dir"; fi
  if [[ -n "$stage" ]]; then rm -rf -- "$stage"; fi
  trap - EXIT
  exit "$status"
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

mkdir -p -- "$stage/raw"; hyb_outputs=()
for ((rep=1; rep<=replicates; rep++)); do
  rep_dir="$stage/raw/replicate-$rep"; mkdir -- "$rep_dir"; rep_id="${run_id}_r${rep}"
  set +e
  (
    cd -- "$rep_dir"
    env LC_ALL=C LANG=C TZ=UTC PYTHONHASHSEED=0 \
      PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0 \
      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      CORES=1 HYB_THREADS=1 \
      "$hyb_bin" detect in="$reads" id="$rep_id" db="$hyb_db" \
        qc=none align=bowtie2 type=mim fold=UNAfold
  ) >"$rep_dir/hyb.stdout.log" 2>"$rep_dir/hyb.stderr.log"
  hyb_status=$?
  set -e
  if [[ "$hyb_status" -ne 0 ]]; then
    cat -- "$rep_dir/hyb.stderr.log" >&2
    echo "Hyb replicate $rep failed with status $hyb_status" >&2
    exit "$hyb_status"
  fi
  hyb_file="$rep_dir/${rep_id}_comp_${hyb_db}_hybrids_ua.hyb"
  if [[ ! -s "$hyb_file" ]]; then
    cat -- "$rep_dir/hyb.stderr.log" >&2
    echo "Hyb did not create non-empty expected output: $hyb_file" >&2
    exit 65
  fi
  hyb_outputs+=("$hyb_file")
done

parser_args=(--hyb "${hyb_outputs[@]}" --sites "$stage/sites.tsv" --targets "$stage/targets.tsv"
  --excluded "$stage/excluded.tsv" --support "$stage/support.tsv" --manifest "$stage/manifest.json"
  --hyb-commit 028ab6371ce793ca5e86f475fce1f2cc6ad3c677 --hyb-db "$hyb_db" --run-id "$run_id")
parser_args+=(--reads "$reads")
if [[ -n "$expressed_mirnas" ]]; then
  parser_args+=(--expression "$expressed_mirnas" --expression-threshold "$expression_threshold")
fi
python3 "$script_dir/consensus_hyb.py" "${parser_args[@]}"
[[ -s "$stage/sites.tsv" && -s "$stage/targets.tsv" && -s "$stage/support.tsv" && -s "$stage/manifest.json" ]] || {
  echo "Parser postcondition failed: missing structured outputs" >&2; exit 65;
}
python3 - "$stage/manifest.json" <<'PY'
import json, sys
with open(sys.argv[1], encoding="utf-8") as handle:
    manifest = json.load(handle)
if manifest.get("status") != "complete" or manifest.get("accepted_rows", 0) < 1:
    raise SystemExit("Parser postcondition failed: no accepted consensus rows")
PY

if [[ -e "$output_dir" ]]; then
  backup="$output_parent/.${output_name}.backup.$$"
  [[ ! -e "$backup" ]] || { echo "Backup path collision: $backup" >&2; exit 73; }
  mv -- "$output_dir" "$backup"
fi
mv -- "$stage" "$output_dir"
stage=
if [[ -n "$backup" ]]; then rm -rf -- "$backup"; backup=; fi
trap - EXIT HUP INT TERM
printf 'Completed consensus Hyb analysis: %s\n' "$output_dir"
