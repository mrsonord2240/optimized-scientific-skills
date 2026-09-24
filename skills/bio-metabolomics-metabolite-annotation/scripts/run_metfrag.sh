#!/usr/bin/env bash
# Rank candidate structures with MetFragCommandLine (LocalCSV candidate list).
# Inputs:  JAR            MetFragCommandLine-<version>.jar (needs Java)
#          PEAKLIST       "mz intensity" per line, no header
#          CANDIDATES     CSV: Identifier,MolecularFormula,MonoisotopicMass,InChI,InChIKey,SMILES,Name
#                         (CSV-quote any field containing a comma)
#          NEUTRAL_MASS   neutral precursor mass of the feature
#          SAMPLE         output name; result is <OUTDIR>/<SAMPLE>.csv, one row per candidate, Score 0-1
#          OUTDIR         optional, default .  (paths are relative to the current directory; on
#                         Windows give Java a path it can read, not an MSYS /c/... path)
# Usage:   bash scripts/run_metfrag.sh JAR PEAKLIST CANDIDATES NEUTRAL_MASS SAMPLE [OUTDIR]
# Example: bash scripts/run_metfrag.sh MetFragCommandLine-2.6.1.jar examples/metfrag/peaklist.txt \
#            examples/metfrag/candidates.csv 192.0270 citrate_test out
# Checked on MetFragCommandLine 2.6.1.
set -euo pipefail
[ $# -ge 5 ] || { sed -n '2,15p' "$0"; exit 2; }
JAR=$1; PEAKLIST=$2; CANDIDATES=$3; NEUTRAL_MASS=$4; SAMPLE=$5; OUTDIR=${6:-.}
mkdir -p "$OUTDIR"

# MetFrag reads "key = value" pairs, one per line; paths are relative to the current directory.
cat > "$OUTDIR/$SAMPLE.params.txt" <<PARAMS
PeakListPath = $PEAKLIST
MetFragDatabaseType = LocalCSV
LocalDatabasePath = $CANDIDATES
NeutralPrecursorMass = $NEUTRAL_MASS
FragmentPeakMatchAbsoluteMassDeviation = 0.01
FragmentPeakMatchRelativeMassDeviation = 10
MetFragCandidateWriter = CSV
SampleName = $SAMPLE
ResultsPath = $OUTDIR
PARAMS
# Swap MetFragDatabaseType to PubChem (and drop LocalDatabasePath) to search PubChem instead.

java -jar "$JAR" "$OUTDIR/$SAMPLE.params.txt"

# MetFrag exits 0 even when it stored nothing (an unquoted InChI comma drops that candidate
# row silently), so check the output row count, not the exit code.
N=$(( $(wc -l < "$OUTDIR/$SAMPLE.csv" 2>/dev/null || echo 1) - 1 ))
[ "$N" -gt 0 ] || { echo "ERROR: 0 candidates stored in $OUTDIR/$SAMPLE.csv (check CSV quoting)" >&2; exit 1; }
echo "candidates in $OUTDIR/$SAMPLE.csv: $N"
