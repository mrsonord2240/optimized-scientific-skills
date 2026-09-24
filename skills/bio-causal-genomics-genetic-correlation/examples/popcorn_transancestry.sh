#!/usr/bin/env bash
set -euo pipefail
# Use ancestry-matched scores and set --num-shared to the true overlap.
popcorn compute -v 1 --bfile pop1 --out pop1_scores
popcorn fit -v 1 --cfile pop1_scores.txt --sfile1 pop1.txt --sfile2 pop2.txt out
