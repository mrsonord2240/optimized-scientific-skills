#!/usr/bin/env python
"""CellPhoneDB v5 statistical (label-permutation) analysis, human data.

Purpose: permutation specificity p-values with multi-subunit complexes scored by the limiting
subunit. Seeded and bit-reproducible only at --threads 1; with threads>1 the pool workers do not
replay the RNG stream (use threads>1 for exploration only). The __main__ guard is required on
Windows because score_interactions=True uses multiprocessing.Pool internally.

Inputs: cellphonedb.zip (cellphonedb-data v5 release), meta TSV (barcode -> cell_type),
        normalized (NOT scaled) counts .h5ad.
Output: CellPhoneDB result tables in --out-dir.
Usage:  python cellphonedb_statistical.py cellphonedb.zip meta.tsv counts_normalized.h5ad \
            [--out-dir cpdb_out] [--iterations 1000] [--threads 1] [--seed 1337]
"""
import argparse

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cpdb_zip')
    p.add_argument('meta_tsv')
    p.add_argument('counts_h5ad')
    p.add_argument('--out-dir', default='cpdb_out')
    p.add_argument('--iterations', type=int, default=1000)   # label-permutation null
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--seed', type=int, default=1337)         # default -1 is unseeded
    args = p.parse_args()

    # threshold=0.1: a gene must be expressed in >=10% of a cluster's cells to count
    # pvalue=0.05 reports per-pair significance
    results = cpdb_statistical_analysis_method.call(
        cpdb_file_path=args.cpdb_zip,
        meta_file_path=args.meta_tsv,
        counts_file_path=args.counts_h5ad,
        counts_data='hgnc_symbol',
        threshold=0.1, iterations=args.iterations, pvalue=0.05, debug_seed=args.seed,
        score_interactions=True, threads=args.threads, output_path=args.out_dir)
    # DEG-driven escape from one-vs-rest: cpdb_degs_analysis_method.call(..., degs_file_path=...)
    return results


if __name__ == '__main__':
    main()
