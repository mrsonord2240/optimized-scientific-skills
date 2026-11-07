"""Run source scCODA composition workflow on annotated h5ad."""
import argparse
import pandas as pd
import scanpy as sc
import tensorflow as tf
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("output_csv")
args = p.parse_args()
adata = sc.read_h5ad(args.input)
obs = adata.obs[["sample", "condition", "cell_type"]]
if obs.isna().any().any():
    p.error("sample, condition and cell_type must not be missing")
meta = obs[["sample", "condition"]].drop_duplicates()
if meta["condition"].nunique() < 2 or meta["sample"].duplicated().any() or (meta.groupby("condition", observed=True).size() < 2).any():
    p.error("Need independent samples and one condition per sample")
counts = pd.crosstab(obs["sample"], obs["cell_type"]).reset_index().merge(meta, on="sample")
data = dat.from_pandas(counts, covariate_columns=["sample", "condition"])
tf.random.set_seed(42)
model = mod.CompositionalAnalysis(data, formula="condition", reference_cell_type="automatic")
result = model.sample_hmc()
result.set_fdr(est_fdr=0.1)
result.summary()
result.credible_effects().to_csv(args.output_csv)
