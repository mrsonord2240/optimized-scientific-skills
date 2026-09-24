"""Deterministic upsetplot example. Requires upsetplot 0.9.0 and pandas >=2.2,<3.

Run: python upset_python.py
Writes three PNG/PDF figure pairs; PDF text uses Type 42 embedding.
"""
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from upsetplot import UpSet, from_contents

if not ("2.2" <= pd.__version__ < "3"):
    raise RuntimeError(f"upsetplot 0.9.0 example requires pandas >=2.2,<3; found {pd.__version__}")

np.random.seed(42)
all_genes = [f"Gene{i}" for i in range(1, 501)]
gene_sets = {
    "Treatment_A": np.random.choice(all_genes, 150, replace=False).tolist(),
    "Treatment_B": np.random.choice(all_genes, 130, replace=False).tolist(),
    "Timepoint_Early": np.random.choice(all_genes, 100, replace=False).tolist(),
    "Timepoint_Late": np.random.choice(all_genes, 180, replace=False).tolist(),
    "Pathway_Response": np.random.choice(all_genes, 90, replace=False).tolist(),
}
core_genes = np.random.choice(all_genes, 25, replace=False).tolist()
for key in ["Treatment_A", "Treatment_B", "Pathway_Response"]:
    gene_sets[key] = sorted(set(gene_sets[key]).union(core_genes))

data = from_contents(gene_sets)  # DataFrame indexed by the Boolean membership columns
data["log2FC"] = np.random.normal(0, 1.5, len(data))
data["pvalue"] = 10 ** np.random.uniform(-5, -0.5, len(data))
data["significant"] = data["pvalue"] < 0.05
df_indexed = data.copy()  # from_contents already returns the membership MultiIndex

def add_count_labels(axes):
    for bar in axes["intersections"].patches:
        height = bar.get_height()
        if height:
            axes["intersections"].annotate(f"{height:g}", (bar.get_x() + bar.get_width() / 2, height),
                                             ha="center", va="bottom", fontsize=7)

def save_upset(upset, filename, size):
    fig = plt.figure(figsize=size)  # Do not create a separate axes under the UpSet layout.
    axes = upset.plot(fig=fig)
    add_count_labels(axes)
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    fig.savefig(Path(filename).with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

basic = UpSet(data, subset_size="count", show_counts=False, sort_by="cardinality",
              sort_categories_by="cardinality", max_subset_rank=20)
save_upset(basic, "upset_basic.png", (12, 8))

customized = UpSet(data, subset_size="count", show_counts=False, sort_by="cardinality",
                   sort_categories_by="cardinality", facecolor="#4DBBD5", element_size=46,
                   max_subset_rank=15)
# Exact Treatment_A-and-Treatment_B-only, rather than every superset containing both.
customized.style_subsets(present=["Treatment_A", "Treatment_B"],
                         absent=["Timepoint_Early", "Timepoint_Late", "Pathway_Response"],
                         facecolor="#E64B35")
save_upset(customized, "upset_customized.png", (14, 8))

with_boxplot = UpSet(df_indexed, subset_size="count", show_counts=False, max_subset_rank=20)
with_boxplot.add_catplot(value="log2FC", kind="box", color="#E64B35")
save_upset(with_boxplot, "upset_with_boxplot.png", (14, 10))

print("Saved upset_basic.png, upset_customized.png, upset_with_boxplot.png")
print("Set sizes:", {name: len(genes) for name, genes in gene_sets.items()})
