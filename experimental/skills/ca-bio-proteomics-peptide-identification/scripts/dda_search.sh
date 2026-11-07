#!/usr/bin/env bash
# End-to-end DDA peptide identification: target-decoy database -> database
# search -> Percolator rescoring -> PSMs at 1% FDR.
#
# Engines: ENGINE=sage (default) or ENGINE=comet. Both write a Percolator .pin,
# so the rescoring and FDR steps below are identical for either.
#
# MZML may hold several space-separated paths (no spaces inside a path): they are
# searched together and rescored in one Percolator run, which is how Percolator
# gets enough PSMs to train (see "Rescore" in references/cli_route.md).
#
# Checked on Sage 0.14.6, Comet 2026.02 rev.2, Percolator 3.09.0,
# OpenMS 3.5.0 DecoyDatabase, on Orbitrap Astral 5-min DDA runs.
set -euo pipefail

MZML="${MZML:?set MZML to one or more centroided mzML paths, space-separated (msconvert --mzML --zlib --filter 'peakPicking vendor msLevel=1-')}"
read -r -a MZMLS <<< "$MZML"
FASTA="${FASTA:?set FASTA to the TARGET-ONLY protein FASTA}"
OUT="${OUT:-./dda_out}"
ENGINE="${ENGINE:-sage}"
THREADS="${THREADS:-8}"
# tool locations -- override to match your install
DECOYDB="${DECOYDB:-DecoyDatabase}"
SAGE="${SAGE:-sage}"
COMET="${COMET:-comet}"
PERCOLATOR="${PERCOLATOR:-percolator}"

mkdir -p "$OUT"

# 1. Target+decoy database. Decoys are made at the PROTEIN level and reversed
#    with the peptide N/C termini held fixed, so decoy peptides obey the same
#    tryptic rules as the targets. Sage generates its own decoys instead;
#    Comet and MS-GF+ need this file.
TDDB="$OUT/target_decoy.fasta"
if [ ! -s "$TDDB" ]; then
  "$DECOYDB" -in "$FASTA" -out "$TDDB" \
    -decoy_string DECOY_ -decoy_string_position prefix \
    -method reverse -enzyme Trypsin -threads "$THREADS"
fi
echo "database: $(grep -c '^>' "$TDDB") entries, $(grep -c '^>DECOY_' "$TDDB") decoys"

case "$ENGINE" in
sage)
  # 2a. Sage generates its own decoys, so it takes the TARGET-ONLY FASTA and
  #     tags decoys 'rev_' (lower case -- see the decoy-prefix trap in SKILL.md).
  # Sage is a native Windows binary in this checked route.  Git Bash does not
  # translate paths embedded in JSON, so convert those paths explicitly; its
  # normal argument conversion still handles the executable invocation.
  SAGE_FASTA="$FASTA"
  SAGE_OUT="$OUT/sage"
  SAGE_MZMLS=("${MZMLS[@]}")
  if command -v cygpath >/dev/null 2>&1; then
    SAGE_FASTA="$(cygpath -m "$FASTA")"
    SAGE_OUT="$(cygpath -m "$OUT/sage")"
    for i in "${!SAGE_MZMLS[@]}"; do SAGE_MZMLS[$i]="$(cygpath -m "${SAGE_MZMLS[$i]}")"; done
  fi
  cat > "$OUT/sage.json" <<JSON
{
  "database": {
    "enzyme": { "missed_cleavages": 2, "min_len": 7, "max_len": 30, "cleave_at": "KR", "restrict": "P" },
    "static_mods": { "C": 57.0215 },
    "variable_mods": { "M": [15.9949] },
    "max_variable_mods": 2,
    "decoy_tag": "rev_",
    "generate_decoys": true,
    "fasta": "$SAGE_FASTA"
  },
  "precursor_tol": { "ppm": [-10, 10] },
  "fragment_tol": { "ppm": [-20, 20] },
  "isotope_errors": [0, 1],
  "deisotope": true,
  "predict_rt": true,
  "report_psms": 1,
  "output_directory": "$SAGE_OUT"
}
JSON
  "$SAGE" "$OUT/sage.json" "${SAGE_MZMLS[@]}" --write-pin  # all runs -> one pin
  PIN="$OUT/sage/results.sage.pin"
  ;;
comet)
  # 2b. Comet searches the concatenated file built in step 1. decoy_search = 0
  #     because the decoys are already in the database; decoy_prefix must match
  #     DecoyDatabase's -decoy_string. Fragment binning 0.02 / offset 0.0 is the
  #     high-res HCD setting; 1.0005 / 0.4 is the ion-trap setting.
  # Like Sage, Comet reads the database and output stem from text rather than
  # shell arguments.  Convert embedded paths under Git Bash so the native
  # executable does not receive a literal /f/... path.
  COMET_TDDB="$TDDB"
  COMET_OUT="$OUT"
  COMET_MZMLS=("${MZMLS[@]}")
  if command -v cygpath >/dev/null 2>&1; then
    COMET_TDDB="$(cygpath -m "$TDDB")"
    COMET_OUT="$(cygpath -m "$OUT")"
    for i in "${!COMET_MZMLS[@]}"; do COMET_MZMLS[$i]="$(cygpath -m "${COMET_MZMLS[$i]}")"; done
  fi
  "$COMET" -p >/dev/null                       # writes comet.params.new
  sed -e "s|^database_name = .*|database_name = $COMET_TDDB|" \
      -e 's|^decoy_search = .*|decoy_search = 0|' \
      -e 's|^decoy_prefix = .*|decoy_prefix = DECOY_|' \
      -e 's|^peptide_mass_tolerance_upper = .*|peptide_mass_tolerance_upper = 10.0|' \
      -e 's|^peptide_mass_tolerance_lower = .*|peptide_mass_tolerance_lower = -10.0|' \
      -e 's|^peptide_mass_units = .*|peptide_mass_units = 2|' \
      -e 's|^isotope_error = .*|isotope_error = 1|' \
      -e 's|^search_enzyme_number = .*|search_enzyme_number = 1|' \
      -e 's|^allowed_missed_cleavage = .*|allowed_missed_cleavage = 2|' \
      -e 's|^fragment_bin_tol = .*|fragment_bin_tol = 0.02|' \
      -e 's|^fragment_bin_offset = .*|fragment_bin_offset = 0.0|' \
      -e 's|^variable_mod01 = .*|variable_mod01 = 15.9949 M 0 2 -1 0 0 0.0|' \
      -e 's|^max_variable_mods_in_peptide = .*|max_variable_mods_in_peptide = 2|' \
      -e 's|^peptide_length_range = .*|peptide_length_range = 7 30|' \
      -e 's|^num_output_lines = .*|num_output_lines = 1|' \
      -e 's|^output_percolatorfile = .*|output_percolatorfile = 1|' \
      -e 's|^output_txtfile = .*|output_txtfile = 1|' \
      comet.params.new > "$OUT/comet.params"
  # One Comet call per run (-N names its output); the pins are then merged with
  # the header kept once, so Percolator trains on all runs together.
  RUN_PINS=()
  n=0
  for f in "${COMET_MZMLS[@]}"; do
    n=$((n + 1))
    "$COMET" -P"$OUT/comet.params" -N"$COMET_OUT/comet_run$n" "$f"
    RUN_PINS+=("$OUT/comet_run$n.pin")
  done
  PIN="$OUT/comet.pin"
  head -n 1 "${RUN_PINS[0]}" > "$PIN"
  for p in "${RUN_PINS[@]}"; do tail -n +2 "$p" >> "$PIN"; done
  ;;
*)
  echo "ENGINE must be sage or comet" >&2
  exit 2
  ;;
esac

# 3. Rescore. Both engines ran ONE concatenated search with one hit per
#    spectrum, so target-decoy competition is the right post-processing, not the
#    mix-max default (-y) that applies to separate target and decoy searches.
#    Pre-flight: a decoy-tag mismatch leaves the pin's Label column with no -1
#    rows, and Percolator would only fail into its log. Stop here, loudly.
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "Label") l = i; next }
            l { n[$l]++ }
            END { printf "pin Label column: %d targets (1), %d decoys (-1)\n", n[1], n[-1]
                  exit !(l && n[1] > 0 && n[-1] > 0) }' "$PIN" \
  || { echo "ERROR: $PIN has no target or no decoy rows -- decoy tag mismatch between the database and the engine's decoy setting (see the decoy-tag table in SKILL.md)" >&2; exit 1; }

if ! "$PERCOLATOR" --post-processing-tdc \
  --results-psms "$OUT/psms.target.tsv" --decoy-results-psms "$OUT/psms.decoy.tsv" \
  --results-peptides "$OUT/peptides.target.tsv" --decoy-results-peptides "$OUT/peptides.decoy.tsv" \
  "$PIN" 2> "$OUT/percolator.log"; then
  echo "ERROR: Percolator failed; last lines of $OUT/percolator.log:" >&2
  tail -n 5 "$OUT/percolator.log" >&2
  exit 1
fi

# 4. The 1% list. Percolator's q-value column is already the list-level FDR;
#    psms.target.tsv holds targets only, so no decoy filtering is needed here.
#    Find q-value BY NAME: Percolator emits a 'filename' column only when the
#    pin has one (Sage does, Comet does not), so its index is not fixed.
cut_1pct() {
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "q-value") q = i; next }
              q && $q <= 0.01' "$1"
}
cut_1pct "$OUT/psms.target.tsv" > "$OUT/psms_1pct.tsv"
cut_1pct "$OUT/peptides.target.tsv" > "$OUT/peptides_1pct.tsv"
echo "PSMs at q <= 0.01:     $(wc -l < "$OUT/psms_1pct.tsv")"
echo "peptides at q <= 0.01: $(wc -l < "$OUT/peptides_1pct.tsv")"
echo "protein grouping and protein-level FDR -> protein-inference"
