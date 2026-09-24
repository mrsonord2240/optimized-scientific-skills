#!/usr/bin/env python
"""MultiVI joint embedding of a mosaic RNA+ATAC design (some cells paired, some RNA-only).

Input:  .h5mu with mod['rna'] and mod['atac'], both over the full cell union, raw counts in
        .X; cells missing a modality carry all-zero rows in that block (MultiVI detects
        presence per cell from whether the block's counts sum to zero).
Output: <out>.h5mu with obsm['X_multivi'] (imputes the missing modality; not measured data)
Usage:  python scripts/multivi_mosaic.py mosaic.h5mu multivi_out [--max-epochs N]
Checked on scvi-tools 1.5.1, mudata 0.4.1.
"""
import argparse

import mudata as md
import scvi

ap = argparse.ArgumentParser()
ap.add_argument('h5mu')
ap.add_argument('out_prefix')
ap.add_argument('--max-epochs', type=int, default=None, help='default: scvi-tools chooses from dataset size')
args = ap.parse_args()

scvi.settings.seed = 0

# No batch_key: the zero-fill mask is not a sequencing batch. Add batch_key only for real
# sequencing batches, as a separate covariate.
mdata = md.read_h5mu(args.h5mu)
scvi.model.MULTIVI.setup_mudata(
    mdata, modalities={'rna_layer': 'rna', 'atac_layer': 'atac'}
)
model = scvi.model.MULTIVI(
    mdata, n_genes=mdata.mod['rna'].n_vars, n_regions=mdata.mod['atac'].n_vars
)
model.train(max_epochs=args.max_epochs)
mdata.obsm['X_multivi'] = model.get_latent_representation()

mdata.write_h5mu(f'{args.out_prefix}.h5mu')
print(f"latent {mdata.obsm['X_multivi'].shape}")
