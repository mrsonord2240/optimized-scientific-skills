#!/bin/bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXAMPLES_DIR="$(cd "${TEST_DIR}/../examples" && pwd)"
ROOT="$(mktemp -d)"
trap 'rm -rf -- "${ROOT}"' EXIT
OUT="${ROOT}/out"
mkdir -p "${OUT}"

for bad in '..' 'ERR/123' '--max-size'; do
    if bash "${EXAMPLES_DIR}/download_single.sh" "${bad}" "${OUT}" >/dev/null 2>&1; then
        echo "download_single accepted unsafe accession: ${bad}" >&2
        exit 1
    fi
    if bash "${EXAMPLES_DIR}/prefetch_large.sh" "${bad}" "${OUT}" >/dev/null 2>&1; then
        echo "prefetch_large accepted unsafe accession: ${bad}" >&2
        exit 1
    fi
    printf '%s\n' "${bad}" > "${ROOT}/accessions.txt"
    bash "${EXAMPLES_DIR}/download_batch.sh" "${ROOT}/accessions.txt" "${OUT}" >/dev/null
    grep -Fx -- "${bad}" "${OUT}/failed.txt" >/dev/null
done

FIXTURE_DIR="${ROOT}/fixture"
mkdir -p "${FIXTURE_DIR}"
cp "${TEST_DIR}/fixtures/download_batch_success.sh" "${FIXTURE_DIR}/download_batch.sh"
cp "${EXAMPLES_DIR}/sra_safety.sh" "${FIXTURE_DIR}/sra_safety.sh"
chmod +x "${FIXTURE_DIR}/download_batch.sh"
# shellcheck source=../examples/ena_fallback.sh
source "${EXAMPLES_DIR}/ena_fallback.sh"

existing="${OUT}/ERR10419835_1.fastq.gz"
printf 'pre-existing\n' > "${existing}"
if run_ena_fallback ERR10419835 "${OUT}" "${FIXTURE_DIR}"; then
    echo 'ENA fallback overwrote an existing output' >&2
    exit 1
fi
test "$(<"${existing}")" = 'pre-existing'
if find "${OUT}" -name '.sra-owned-stage' -print -quit | grep -q .; then
    echo 'Owned ENA stage was not cleaned after existing-file rejection' >&2
    exit 1
fi
rm -f -- "${existing}"

ln -s "${ROOT}/missing-target" "${existing}"
if run_ena_fallback ERR10419835 "${OUT}" "${FIXTURE_DIR}"; then
    echo 'ENA fallback replaced a dangling symlink' >&2
    exit 1
fi
test -L "${existing}"
if find "${OUT}" -name '.sra-owned-stage' -print -quit | grep -q .; then
    echo 'Owned ENA stage was not cleaned after dangling-symlink rejection' >&2
    exit 1
fi

echo 'SRA safety regressions passed'
