#!/bin/bash
# Reference: DIA-NN README 1.9.2 / 2.0 / 2.6.1 | Verify API if version differs
# DIA-NN 2.x predicted-library analysis: generate from FASTA with NO raw files, then search raw data.
# Filtering the parquet report is done afterwards (points listed at the end; code in SKILL.md).

set -e

FASTA="uniprot_human_reviewed.fasta"
OUTPUT_DIR="diann_out"
PREDICTED_LIBRARY_TEMPLATE="$OUTPUT_DIR/human"
PREDICTED_LIBRARY="$OUTPUT_DIR/human.predicted.speclib"
THREADS=8

# QVALUE 0.01 = 1% precursor FDR (run context); DIA-NN's own default is 0.05.
# MISSED_CLEAVAGES 1 = trypsin standard; higher expands the search/FDR burden.
# MASS_ACC / MASS_ACC_MS1 = fixed per instrument (ppm): timsTOF 15/15, Orbitrap Astral 10/4, TripleTOF 20/20.
#   0 = auto, optimised on the FIRST run and reused for all others, so results depend on run order.
# --reanalyse = MBR: second pass with the empirical library built in the first pass.
# --fixed-mod UniMod:4,57.021464,C = fixed Cys carbamidomethylation (explicit form of the --unimod4 preset).
# --var-mod UniMod:35,15.994915,M = variable Met oxidation.
QVALUE=0.01
MISSED_CLEAVAGES=1
MASS_ACC=15
MASS_ACC_MS1=15

command -v diann >/dev/null || { echo "diann not found on PATH" >&2; exit 1; }
[ -f "$FASTA" ] || { echo "FASTA not found: $FASTA" >&2; exit 1; }
shopt -s nullglob
MZML=(*.mzML)
[ ${#MZML[@]} -gt 0 ] || { echo "No .mzML files in $(pwd)" >&2; exit 1; }
MZML_ARGS=()
for f in "${MZML[@]}"; do
    MZML_ARGS+=(--f "$f")    # array keeps file names with spaces intact
done

mkdir -p "$OUTPUT_DIR"

echo "Running DIA-NN predicted-library analysis..."
# Stage 1: no --f inputs. DIA-NN 2.x requires prediction to be separate from raw-data analysis.
diann \
    --fasta "$FASTA" --fasta-search \
    --gen-spec-lib --predictor \
    --out-lib "$PREDICTED_LIBRARY_TEMPLATE" \
    --cut "K*,R*" --missed-cleavages $MISSED_CLEAVAGES \
    --min-pep-len 7 --max-pep-len 30 \
    --fixed-mod UniMod:4,57.021464,C --var-mods 1 --var-mod UniMod:35,15.994915,M \
    --threads $THREADS

[ -f "$PREDICTED_LIBRARY" ] || { echo "Predicted library was not written: $PREDICTED_LIBRARY" >&2; exit 1; }

# Stage 2: search the raw files using the predicted library. Do not reactivate FASTA digest or prediction.
diann \
    "${MZML_ARGS[@]}" \
    --lib "$PREDICTED_LIBRARY" --fasta "$FASTA" \
    --out "$OUTPUT_DIR/report.parquet" \
    --out-lib "$OUTPUT_DIR/report-lib.parquet" \
    --qvalue $QVALUE --matrices \
    --mass-acc $MASS_ACC --mass-acc-ms1 $MASS_ACC_MS1 \
    --reanalyse --smart-profiling \
    --threads $THREADS

echo "Done. Main report (1.9+ default): $OUTPUT_DIR/report.parquet"
echo "Matrices: $OUTPUT_DIR/report.pg_matrix.tsv (verify *_matrix dotting vs installed version)"

# Filter the parquet report. Points to note for downstream code:
#   - read report.parquet, NOT report.tsv (1.9+ default).
#   - filter BOTH levels: Q.Value (precursor) AND PG.Q.Value (protein-group).
#   - for a cross-run matrix add Global.Q.Value <= 0.01 and Global.PG.Q.Value <= 0.01 (per-run FDR inflates
#     across runs); DIA-NN 1.9.x with MBR: Lib.Q.Value / Lib.PG.Q.Value instead of the Global.* pair.
#   - do not infer FDR or a count direction from report vs matrix; inspect report.log.txt for the installed version's filters.
#   - convert DIA-NN's 0 (not-quantified) to NA before log2 / normalization.
