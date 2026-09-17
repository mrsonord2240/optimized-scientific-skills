#!/usr/bin/env bash
# End-to-end DDA peptide identification: target-decoy database -> database
# search -> Percolator rescoring -> PSMs at 1% FDR.
#
# Engines: ENGINE=sage (default) or ENGINE=comet. Both write a Percolator .pin,
# so the rescoring and FDR steps below are identical for either.
#
# Checked on Sage 0.14.6, Comet 2026.02 rev.2, Percolator 3.09.0,
# OpenMS 3.5.0 DecoyDatabase, on one Orbitrap Astral 5-min DDA run.
set -euo pipefail

MZML="${MZML:?set MZML to the centroided mzML (msconvert --mzML --zlib --filter 'peakPicking vendor msLevel=1-')}"
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
#    tryptic rules as the targets. Sage and MSFragger generate their own decoys
#    instead; Comet, MS-GF+ and X!Tandem need this file.
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
  cat > "$OUT/sage.json" <<JSON
{
  "database": {
    "enzyme": { "missed_cleavages": 2, "min_len": 7, "max_len": 30, "cleave_at": "KR", "restrict": "P" },
    "static_mods": { "C": 57.0215 },
    "variable_mods": { "M": [15.9949] },
    "max_variable_mods": 2,
    "decoy_tag": "rev_",
    "generate_decoys": true,
    "fasta": "$FASTA"
  },
  "precursor_tol": { "ppm": [-10, 10] },
  "fragment_tol": { "ppm": [-20, 20] },
  "isotope_errors": [0, 1],
  "deisotope": true,
  "predict_rt": true,
  "report_psms": 1,
  "output_directory": "$OUT/sage"
}
JSON
  "$SAGE" "$OUT/sage.json" "$MZML" --write-pin
  PIN="$OUT/sage/results.sage.pin"
  ;;
comet)
  # 2b. Comet searches the concatenated file built in step 1. decoy_search = 0
  #     because the decoys are already in the database; decoy_prefix must match
  #     DecoyDatabase's -decoy_string. Fragment binning 0.02 / offset 0.0 is the
  #     high-res HCD setting; 1.0005 / 0.4 is the ion-trap setting.
  "$COMET" -p >/dev/null                       # writes comet.params.new
  sed -e "s|^database_name = .*|database_name = $TDDB|" \
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
  "$COMET" -P"$OUT/comet.params" -N"$OUT/comet" "$MZML"
  PIN="$OUT/comet.pin"
  ;;
*)
  echo "ENGINE must be sage or comet" >&2
  exit 2
  ;;
esac

# 3. Rescore. Both engines ran ONE concatenated search with one hit per
#    spectrum, so target-decoy competition is the right post-processing, not the
#    mix-max default (-y) that applies to separate target and decoy searches.
"$PERCOLATOR" --post-processing-tdc \
  --results-psms "$OUT/psms.target.tsv" --decoy-results-psms "$OUT/psms.decoy.tsv" \
  --results-peptides "$OUT/peptides.target.tsv" --decoy-results-peptides "$OUT/peptides.decoy.tsv" \
  "$PIN" 2> "$OUT/percolator.log"

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
