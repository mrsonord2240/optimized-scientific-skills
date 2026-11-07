#!/usr/bin/env python3
'''
Single-cell splicing analysis using BRIE2.
Estimates PSI values with uncertainty quantification for sparse scRNA-seq data.
'''
# Reference: checked on brie 2.3.0 (GitHub install), anndata 0.12, numpy 2.x, pandas 2.3, scanpy 1.11 | Verify API if version differs

import shutil
import subprocess
import urllib.request

ANNOTATION_URL = 'https://sourceforge.net/projects/brie-rna/files/annotation/{path}/SE.{tier}.gff3/download'
ANNOTATION_PATH = {'human': 'human/gencode.v25', 'mouse': 'mouse/gencode.vM12'}


def fetch_splicing_events(species, output_gff, tier='most'):
    '''
    Download BRIE's precomputed exon-skipping (SE) event annotation.

    species: 'human' (GENCODE v25) or 'mouse' (GENCODE vM12); align reads to the same or a
        close genome version. tier: 'gold' (strictest) or 'most' (more events).
    The briekit-event generator is not used: its entry point crashes on install
    (ModuleNotFoundError: parseTables).
    '''
    url = ANNOTATION_URL.format(path=ANNOTATION_PATH[species], tier=tier)
    with urllib.request.urlopen(url, timeout=120) as resp, open(output_gff, 'wb') as out:
        shutil.copyfileobj(resp, out)

    print(f'Splicing events written to: {output_gff}')


def count_splicing_reads(bam_file, events_gff, output_dir, barcode_file, n_proc=15):
    '''
    Count splice-junction reads per cell with brie-count (10X droplet mode).

    Args:
        bam_file: Possorted BAM from Cell Ranger
        events_gff: Splicing-event GFF3 from fetch_splicing_events
        output_dir: Output directory (writes brie_count.h5ad)
        barcode_file: Filtered barcodes (barcodes.tsv[.gz])
        n_proc: Parallel processes
    '''
    subprocess.run([
        'brie-count', '-a', events_gff, '-s', bam_file,
        '-b', barcode_file, '-o', output_dir, '-p', str(n_proc)
    ], check=True)

    print(f'BRIE counts written to: {output_dir}/brie_count.h5ad')


def run_brie2_inference(count_h5ad, cell_features, output_h5ad):
    '''
    Run BRIE2 quantification (brie-quant) for PSI + cell-feature effects.

    Args:
        count_h5ad: brie_count.h5ad from brie-count
        cell_features: TSV of cell covariates for the design matrix
        output_h5ad: Output AnnData path

    Returns AnnData with PSI estimates in layers ('Psi').
    '''
    import brie

    # Probabilistic PSI estimation handles sparse data
    subprocess.run([
        'brie-quant', '-i', count_h5ad, '-c', cell_features,
        '-o', output_h5ad, '--interceptMode', 'gene', '--LRTindex', 'All'
    ], check=True)

    return brie.read_h5ad(output_h5ad)


def find_variable_splicing(adata_splice, cell_type_col='cell_type', min_cells=50):
    '''
    Identify splicing events variable across cell types.

    Args:
        adata_splice: AnnData from BRIE2 with PSI in layers
        cell_type_col: Column in obs with cell type labels
        min_cells: Minimum cells per type for reliable estimates
    '''
    import numpy as np
    import pandas as pd

    psi = adata_splice.layers['Psi']
    cell_types = adata_splice.obs[cell_type_col].unique()

    mean_psi = pd.DataFrame(index=adata_splice.var_names)

    for ct in cell_types:
        mask = adata_splice.obs[cell_type_col] == ct
        if mask.sum() >= min_cells:
            mean_psi[ct] = np.nanmean(psi[mask, :], axis=0)

    # Calculate variance across cell types
    psi_variance = mean_psi.var(axis=1)
    variable_events = psi_variance.sort_values(ascending=False)

    return variable_events, mean_psi


def pseudobulk_by_celltype(adata, cell_type_col='cell_type', sample_col='sample',
                           min_cells=10, min_replicates=3, layer=None):
    '''
    Create pseudobulk samples by summing counts per (sample, cell type).
    Useful for running bulk splicing tools on scRNA-seq.

    One pooled column per cell type is n = 1 per group and gives false positives, so
    every cell type needs >= min_replicates samples (donors/batches) with >= min_cells cells.
    Counts come from adata.X, or adata.layers[layer] (e.g. 'isoform2' from brie-count).
    '''
    import numpy as np
    import pandas as pd
    from scipy import sparse

    X = adata.layers[layer] if layer else adata.X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)

    obs = adata.obs
    if obs[[cell_type_col, sample_col]].isna().any().any():
        raise ValueError(f'NaN in obs[{cell_type_col!r}] or obs[{sample_col!r}]: those cells would be dropped')

    pseudobulk = {}
    n_rep = {}

    for (ct, sample), idx in obs.groupby([cell_type_col, sample_col], observed=True).indices.items():
        if len(idx) >= min_cells:
            pseudobulk[f'{ct}__{sample}'] = X[idx].sum(axis=0)
            n_rep[ct] = n_rep.get(ct, 0) + 1

    too_few = {ct: n for ct, n in n_rep.items() if n < min_replicates}
    if too_few or set(n_rep) != set(obs[cell_type_col].unique()):
        raise ValueError(f'need >= {min_replicates} samples with >= {min_cells} cells per cell type; got {n_rep}')

    return pd.DataFrame(pseudobulk, index=adata.var_names)


def differential_splicing_per_cell(adata_splice, group1_cells, group2_cells):
    '''
    Compare per-cell Psi between two cell populations with a Mann-Whitney test.

    Cells are treated as replicates, so cells from the same donor are pseudo-replicates;
    this is not a pseudobulk test (use pseudobulk_by_celltype + a bulk tool for that).

    Args:
        adata_splice: AnnData with PSI estimates
        group1_cells: Boolean mask or list of cell barcodes for group 1
        group2_cells: Boolean mask or list of cell barcodes for group 2
    '''
    import numpy as np
    from scipy import stats

    psi = adata_splice.layers['Psi']

    if isinstance(group1_cells, list):
        mask1 = adata_splice.obs_names.isin(group1_cells)
        mask2 = adata_splice.obs_names.isin(group2_cells)
    else:
        mask1 = group1_cells
        mask2 = group2_cells

    results = []

    for i, event in enumerate(adata_splice.var_names):
        psi1 = psi[mask1, i]
        psi2 = psi[mask2, i]

        # Remove NaN values
        psi1 = psi1[~np.isnan(psi1)]
        psi2 = psi2[~np.isnan(psi2)]

        if len(psi1) >= 20 and len(psi2) >= 20:
            # Mann-Whitney U test
            stat, pval = stats.mannwhitneyu(psi1, psi2, alternative='two-sided')
            delta_psi = np.mean(psi1) - np.mean(psi2)

            results.append({
                'event': event,
                'mean_psi_group1': np.mean(psi1),
                'mean_psi_group2': np.mean(psi2),
                'delta_psi': delta_psi,
                'pvalue': pval,
                'n_cells_group1': len(psi1),
                'n_cells_group2': len(psi2)
            })

    import pandas as pd
    results_df = pd.DataFrame(results)

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    if len(results_df) > 0:
        _, results_df['fdr'], _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
        results_df = results_df.sort_values('fdr')

    return results_df


if __name__ == '__main__':
    # Example workflow
    print('BRIE2 Single-Cell Splicing Analysis')
    print('=' * 40)

    # Step 1: Get events
    # fetch_splicing_events('mouse', 'splicing_events.gff3')

    # Step 2: Count reads
    # count_splicing_reads(
    #     'possorted_genome_bam.bam',
    #     'splicing_events.gff3',
    #     'brie_counts/',
    #     'filtered_barcodes.tsv'
    # )

    # Step 3: Run BRIE2 inference
    # adata_splice = run_brie2_inference(
    #     'brie_counts/brie_count.h5ad', 'cell_features.tsv', 'brie_quant.h5ad'
    # )

    # Step 4: Find variable events
    # variable_events, mean_psi = find_variable_splicing(adata_splice)
    # print('Top variable splicing events:')
    # print(variable_events.head(20))

    print('Provide input files to run single-cell splicing analysis')
