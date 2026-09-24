"""Create a tiny deterministic AnnData smoke-test fixture for this skill.

It checks file I/O, required obs fields, gene-symbol handling, and the Windows-safe
runner before users spend time on a full dataset. It is not biological validation.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def make_fixture(output_path):
    rng = np.random.default_rng(12345)
    genes = ["NT5E", "PTGES", "PTGER2", "PTGER4", "SLC16A1", "GAPDH"]
    groups = ["Tumor"] * 12 + ["TCell"] * 12
    matrix = rng.poisson(1.0, size=(24, len(genes))).astype(float)
    matrix[:12, genes.index("NT5E")] += 6
    matrix[:12, genes.index("PTGES")] += 6
    matrix[12:, genes.index("PTGER2")] += 6
    matrix[12:, genes.index("PTGER4")] += 6
    adata = ad.AnnData(
        X=np.log1p(matrix),
        obs=pd.DataFrame({"cell_type": groups}, index=[f"cell_{i}" for i in range(24)]),
        var=pd.DataFrame(index=genes),
    )
    adata.write_h5ad(output_path)
    return adata


if __name__ == "__main__":
    output = Path(__file__).with_name("synthetic_metabolite_communication.h5ad")
    adata = make_fixture(output)
    print(f"Wrote {output.name}: {adata.n_obs} cells x {adata.n_vars} genes")
