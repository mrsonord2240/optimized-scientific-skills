#!/bin/bash
set -euo pipefail

accessions="$1"
out_dir="$2"
accession="$(head -n 1 "${accessions}")"
printf 'verified fixture mate 1\n' > "${out_dir}/${accession}_1.fastq.gz"
printf 'verified fixture mate 2\n' > "${out_dir}/${accession}_2.fastq.gz"
: > "${out_dir}/failed.txt"
