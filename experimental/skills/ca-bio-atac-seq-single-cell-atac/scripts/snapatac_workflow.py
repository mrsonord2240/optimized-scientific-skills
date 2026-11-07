#!/usr/bin/env python3
"""Source SnapATAC2 workflow on indexed hg38 fragments."""
import argparse
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument('fragments')
    p.add_argument('output')
    a = p.parse_args()
    if not Path(a.fragments).is_file():
        p.error('fragment file does not exist')
    if Path(a.output).exists():
        p.error('output exists; choose a fresh path')
    import snapatac2 as snap
    data = snap.pp.import_fragments(fragment_file=a.fragments, chrom_sizes=snap.genome.hg38,
                                    file=a.output, sorted_by_barcode=False)
    snap.metrics.tsse(data, gene_anno=snap.genome.hg38)
    snap.metrics.frag_size_distr(data)
    snap.pp.filter_cells(data, min_counts=1000, min_tsse=4)
    snap.pp.add_tile_matrix(data, bin_size=500)
    snap.pp.select_features(data, n_features=250000)
    snap.tl.spectral(data)
    snap.pp.knn(data)
    snap.tl.umap(data)
    snap.tl.leiden(data)
    snap.tl.macs3(data, groupby='leiden', n_jobs=1)

if __name__ == '__main__':
    main()
