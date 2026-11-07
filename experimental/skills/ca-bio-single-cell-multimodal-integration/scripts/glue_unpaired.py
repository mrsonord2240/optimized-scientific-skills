#!/usr/bin/env python3
"""GLUE alignment of unpaired scRNA and scATAC through a peak-near-gene guidance graph.

Usage:  python glue_unpaired.py rna.h5ad atac.h5ad out_prefix
Input:  rna.h5ad: RAW counts in X (NB likelihood), highly variable genes flagged, X_pca in obsm, gene coordinates
        in var (chrom, chromStart, chromEnd, strand). atac.h5ad: counts in X (ZINB), X_lsi in obsm, peak names
        chr:start-end in var. RNA and ATAC coordinates must share one genome build or the guidance graph is
        near-empty and the alignment is garbage.
Output: out_prefix_rna.h5ad and out_prefix_atac.h5ad with obsm['X_glue'] (aligned embeddings).
After running, verify cell-type structure is preserved, not only modality overlap: adversarial alignment can
over-mix distinct populations.
Needs:  scglue. No Windows install (pybedtools -> pysam has no wheel): use WSL, Linux or macOS.
Checked: against scglue's documented API only; not run, not yet executed after the recut.
"""
import sys

import anndata as ad
import scglue

if len(sys.argv) != 4:
    sys.exit(__doc__)
rna = ad.read_h5ad(sys.argv[1])
atac = ad.read_h5ad(sys.argv[2])
out = sys.argv[3]

scglue.models.configure_dataset(rna, 'NB', use_highly_variable=True, use_rep='X_pca')
scglue.models.configure_dataset(atac, 'ZINB', use_highly_variable=True, use_rep='X_lsi')
graph = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
print('guidance graph edges:', graph.number_of_edges())  # near 0 means a genome-build or coordinate mismatch
# GLUE trains a VAE; the seed is pinned (scglue's documented default is already random_seed=0).
glue = scglue.models.fit_SCGLUE({'rna': rna, 'atac': atac}, graph, init_kws={'random_seed': 0})
rna.obsm['X_glue'] = glue.encode_data('rna', rna)
atac.obsm['X_glue'] = glue.encode_data('atac', atac)
rna.write(f'{out}_rna.h5ad')
atac.write(f'{out}_atac.h5ad')
print('wrote', f'{out}_rna.h5ad', f'{out}_atac.h5ad')
