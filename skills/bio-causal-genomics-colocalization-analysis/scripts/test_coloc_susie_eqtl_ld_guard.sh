#!/usr/bin/env bash
# test_coloc_susie_eqtl_ld_guard.sh -- execute the eQTL LD guard regression on Windows.
# Inputs: R wrapper, test_coloc_susie_eqtl_ld_guard.R, coloc_susie.R.
# Usage:  bash scripts/test_coloc_susie_eqtl_ld_guard.sh <r.sh> <test.R> <coloc_susie.R>

set -euo pipefail
if [[ "$#" -ne 3 ]]; then
  echo "usage: $0 <r.sh> <test.R> <coloc_susie.R>" >&2
  exit 2
fi

set +e
output=$("$1" "$2" "$3" 2>&1)
set -e
printf '%s\n' "$output"

# The R regression intentionally invokes a child Rscript that stops on the bad
# eQTL LD panel. This runtime propagates that child status to the parent, so the
# PASS marker emitted only after all R assertions is the test result.
grep -Fqx 'PASS: eQTL LD mismatch rejected (lambda=1.000000) before coloc output' <<< "$output"
