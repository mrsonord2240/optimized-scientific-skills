#!/bin/bash
# Smoke test: confirms an LDSC install actually works before pointing it at real GWAS data.
# Exercises --h2, --rg, and --h2-cts (Finucane 2018 cell-type prioritization) against the tiny
# simulated fixtures CBIIT/ldsc ships in its own test/ directory (1000-SNP Z-score panels and
# chromosome-split LD scores; not real GWAS data). They are read from your clone, not bundled here.
#
# Prerequisite: git clone https://github.com/CBIIT/ldsc.git, Python 3.9+ env per SKILL.md's
# "Tool Install Notes", and the one-line ldscore/sumstats.py patch (.loc[:,1:] -> .iloc[:,1:])
# for the --h2-cts step below.
#
# Usage:
#   LDSC_DIR=/path/to/ldsc bash smoke_test_ldsc.sh

set -euo pipefail

LDSC_DIR=$(cd "${LDSC_DIR:?Set LDSC_DIR to your CBIIT/ldsc clone, e.g. LDSC_DIR=./ldsc bash smoke_test_ldsc.sh}" && pwd)
OUT=$(mktemp -d)

# Stage the fixtures. CBIIT ships one unsplit regression-weight file; --h2-cts reads weights split
# by chromosome, so the same file is copied per chromosome -- a smoke-test convenience only.
cp -r "${LDSC_DIR}/test/simulate_test/ldscore" "${LDSC_DIR}/test/simulate_test/sumstats" "${OUT}/"
cp "${OUT}/ldscore/w.l2.ldscore" "${OUT}/ldscore/w1.l2.ldscore"
cp "${OUT}/ldscore/w.l2.ldscore" "${OUT}/ldscore/w2.l2.ldscore"
gzip "${OUT}"/ldscore/*.l2.ldscore  # ldsc.py opens LD score files only as .gz
# Two two-category LD score sets standing in for "cell types"; there is no correct tissue here.
printf 'CellTypeA\tldscore/twold_firstfile\nCellTypeB\tldscore/twold_secondfile\n' > "${OUT}/test.ldcts"

echo "=== munge_sumstats.py ==="
python "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${LDSC_DIR}/test/munge_test/sumstats" \
    --merge-alleles "${LDSC_DIR}/test/munge_test/merge_alleles" \
    --N 6702 --signed-sumstats OR,1 \
    --out "${OUT}/munge"
# Expect: "1 SNPs remain" and a printed Mean chi^2 / Lambda GC -- confirms munging runs.

echo "=== --h2 (total heritability) ==="
python "${LDSC_DIR}/ldsc.py" --h2 "${OUT}/sumstats/0" \
    --ref-ld "${OUT}/ldscore/oneld_onefile1" \
    --w-ld "${OUT}/ldscore/w1" \
    --out "${OUT}/h2"
grep "Total Observed scale h2" "${OUT}/h2.log"
# Expect: Total Observed scale h2 ~ 0.38 (0.04), Intercept ~2.09, Ratio ~0.12.

echo "=== --rg (cross-trait genetic correlation) ==="
python "${LDSC_DIR}/ldsc.py" --rg "${OUT}/sumstats/0,${OUT}/sumstats/1" \
    --ref-ld "${OUT}/ldscore/twold_onefile1" \
    --w-ld "${OUT}/ldscore/w1" \
    --out "${OUT}/rg"
grep "Genetic Correlation" "${OUT}/rg.log"
# Expect: Genetic Correlation ~ 0.11 (0.08), P ~ 0.15 (these two simulated traits are ~uncorrelated).

echo "=== --h2-cts (Finucane 2018 cell-type prioritization) ==="
echo "(requires the ldscore/sumstats.py .loc[:,1:] -> .iloc[:,1:] patch -- see SKILL.md)"
# test.ldcts' rows are relative paths, so run from the staging directory.
( cd "${OUT}" && python "${LDSC_DIR}/ldsc.py" --h2-cts "sumstats/0" \
    --ref-ld-chr "ldscore/oneld_onefile" \
    --ref-ld-chr-cts "test.ldcts" \
    --w-ld-chr "ldscore/w" \
    --out "${OUT}/cts" )
cat "${OUT}/cts.cell_type_results.txt"
# Expect a two-row table (CellTypeA, CellTypeB) with Coefficient / SE / P-value: one row well under
# a Bonferroni threshold of 0.05/2 = 0.025 (P ~ 3.5e-4), the other not (P ~ 0.99).

echo "All three LDSC entry points ran. Output in ${OUT}/"
