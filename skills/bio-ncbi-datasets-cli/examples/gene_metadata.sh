#!/bin/bash
# Reference: NCBI Datasets CLI 18.37.0 (checked 2026-09-19) | Verify API if version differs
# Cross-species gene metadata via Datasets summary + dataformat TSV.
#
# `--taxon` on `summary gene symbol` is single-species only -- it does NOT accept a
# clade like "Mammalia" (errors: "gene requires an at-or-below-species-level taxon").
# The only cross-species mechanism on this subcommand is `--ortholog <taxon|all>`,
# which returns NCBI's ortholog set (one representative gene per species, limited to
# vertebrates and insects) filtered to the given clade -- so ORTHOLOG_TAXON below can
# be any rank ("Mammalia", "Insecta", a species name, ...) or the literal "all".

set -euo pipefail

SYMBOL="${1:-BRCA1}"
ORTHOLOG_TAXON="${2:-Mammalia}"

echo "=== ${SYMBOL} gene records across ${ORTHOLOG_TAXON} (NCBI ortholog set, one rep per species) ==="
datasets summary gene symbol "${SYMBOL}" \
    --ortholog "${ORTHOLOG_TAXON}" \
    --as-json-lines \
  | dataformat tsv gene \
        --fields gene-id,symbol,tax-name,description,chromosomes \
  > "${SYMBOL}_${ORTHOLOG_TAXON// /_}.tsv"

head "${SYMBOL}_${ORTHOLOG_TAXON// /_}.tsv" | column -t -s $'\t'

echo
echo "Note: NCBI's ortholog set is one representative gene per species, limited to"
echo "vertebrates and insects. For tree-reconciled orthology with co-orthologs and"
echo "1:many calls, use ortholog-inference (Compara, OMA) instead."

echo
echo "=== Single-species lookup (the other supported --taxon usage) ==="
datasets summary gene symbol "${SYMBOL}" --taxon human --as-json-lines \
  | dataformat tsv gene \
        --fields gene-id,symbol,tax-name,description \
  > "${SYMBOL}_human.tsv"

head "${SYMBOL}_human.tsv" | column -t -s $'\t'

echo
echo "=== Inspect raw JSON-lines for fields available ==="
echo "datasets summary gene symbol ${SYMBOL} --taxon human --as-json-lines | jq -s 'first | keys'"
echo "Use field names from there, or 'dataformat tsv gene --help', in dataformat --fields"
