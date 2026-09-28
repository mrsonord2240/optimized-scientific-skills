"""Compare PCA, UMAP, t-SNE, and PHATE on a raw-count AnnData file.

Input contract: at least 100 cells and 2,000 genes; ``adata.X`` contains raw,
non-negative integer counts; ``adata.obs['condition']`` is present. Optional
``adata.obs['pseudotime']`` colors PHATE.

Usage:
    python embedding_phd.py raw_counts.h5ad --output-dir figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openTSNE
import phate
import scanpy as sc
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_h5ad", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("figures"))
    return parser.parse_args()


def validate_input(adata: sc.AnnData) -> None:
    if adata.n_obs < 100 or adata.n_vars < 2_000:
        raise ValueError("input requires at least 100 cells and 2,000 genes")
    if "condition" not in adata.obs:
        raise ValueError("adata.obs['condition'] is required")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if np.any(values < 0) or not np.allclose(values, np.round(values)):
        raise ValueError("adata.X must contain raw, non-negative integer counts")


def neighbor_retention(high_dim: np.ndarray, embedding: np.ndarray, k: int = 15) -> float:
    """Return mean overlap between high-dimensional and 2D k-neighbor sets."""
    if len(high_dim) <= k:
        raise ValueError(f"neighbor retention requires more than k={k} rows")

    def indices(values: np.ndarray) -> np.ndarray:
        model = NearestNeighbors(n_neighbors=k + 1).fit(values)
        return model.kneighbors(values, return_distance=False)[:, 1:]

    high_neighbors = indices(high_dim)
    low_neighbors = indices(embedding)
    return float(
        np.mean(
            [len(set(high) & set(low)) / k for high, low in zip(high_neighbors, low_neighbors)]
        )
    )


def categorical_colors(n_categories: int) -> np.ndarray:
    """Return one distinct tab20 RGBA value per category, up to twenty."""
    if not 1 <= n_categories <= 20:
        raise ValueError(
            "categorical color encoding supports 1-20 categories; "
            "facet the plot or supply a documented alternative encoding above 20"
        )
    return np.asarray(plt.colormaps["tab20"](np.arange(n_categories)))


def annotate_loading_labels(
    ax: plt.Axes, endpoints: np.ndarray, labels: list[str]
) -> list[plt.Annotation]:
    """Place loading labels with deterministic, rendered-box collision checks."""
    if len(endpoints) != len(labels):
        raise ValueError("loading endpoints and labels must have equal length")

    # Search close to each arrow tip first, alternating below and above. Wider
    # columns provide a deterministic escape hatch for unusually long labels.
    offsets = [
        (x_offset, y_offset)
        for x_offset in (6, 30, 54)
        for y_offset in (0, -12, 12, -24, 24, -36, 36, -48, 48)
    ]
    accepted_boxes = []
    annotations = []
    for endpoint, label in zip(endpoints, labels):
        for offset in offsets:
            annotation = ax.annotate(
                label,
                xy=tuple(endpoint),
                xytext=offset,
                textcoords="offset points",
                fontsize=6,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.7,
                    "edgecolor": "none",
                    "pad": 0.2,
                },
            )
            ax.figure.canvas.draw()
            box = annotation.get_window_extent(ax.figure.canvas.get_renderer())
            if not any(box.overlaps(accepted) for accepted in accepted_boxes):
                annotations.append(annotation)
                accepted_boxes.append(box)
                break
            annotation.remove()
        else:
            raise RuntimeError("could not place PCA loading labels without overlap")
    return annotations


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sc.set_figure_params(dpi_save=300, figsize=(4, 4), frameon=False)

    adata = sc.read_h5ad(args.input_h5ad)
    validate_input(adata)

    # seurat_v3 expects raw counts, so select HVGs before normalization/log1p.
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2_000, flavor="seurat_v3", subset=True
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, random_state=42)

    scores = adata.obsm["X_pca"]
    variance = adata.uns["pca"]["variance_ratio"]
    conditions = adata.obs["condition"].astype("category")
    categories = list(conditions.cat.categories)
    condition_colors = categorical_colors(len(categories))

    # PCA with truthful categorical legend and the five strongest PC1/PC2 loadings.
    fig, ax = plt.subplots(figsize=(5, 4))
    for index, category in enumerate(categories):
        selected = np.asarray(conditions == category)
        ax.scatter(
            scores[selected, 0], scores[selected, 1], color=condition_colors[index],
            label=str(category), alpha=0.7, s=10, rasterized=True,
        )
    loadings = adata.varm["PCs"][:, :2]
    strongest = np.argsort(np.linalg.norm(loadings, axis=1))[-5:]
    arrow_scale = 0.12 * min(np.ptp(scores[:, 0]), np.ptp(scores[:, 1]))
    loading_scale = arrow_scale / max(np.linalg.norm(loadings[strongest], axis=1))
    endpoints = []
    loading_labels = []
    for gene_index in strongest:
        dx, dy = loadings[gene_index] * loading_scale
        ax.arrow(0, 0, dx, dy, color="black", width=0.002 * arrow_scale)
        endpoints.append((dx, dy))
        loading_labels.append(str(adata.var_names[gene_index]))
    annotate_loading_labels(ax, np.asarray(endpoints), loading_labels)
    ax.set_xlabel(f"PC1 ({variance[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance[1] * 100:.1f}%)")
    ax.set_title("PCA")
    ax.legend(title="condition", bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, args.output_dir / "pca.pdf")

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(range(1, len(variance) + 1), variance, "o-", markersize=3)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained")
    save_figure(fig, args.output_dir / "scree.pdf")

    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.3, random_state=42)
    sc.tl.leiden(
        adata, resolution=0.5, random_state=42,
        flavor="igraph", n_iterations=2, directed=False,
    )
    umap_ax = sc.pl.umap(
        adata, color="leiden", palette="tab20", legend_loc="on data",
        legend_fontsize=7, show=False,
    )
    save_figure(umap_ax.figure, args.output_dir / "umap_clusters.pdf")

    # openTSNE 1.0 already defaults to PCA initialization and auto=n/12 learning
    # rate; spell them out so the recorded analysis contract is explicit.
    tsne_embedding = openTSNE.TSNE(
        perplexity=30,
        n_iter=750,
        initialization="pca",
        learning_rate="auto",
        n_jobs=-1,
        random_state=42,
    ).fit(scores[:, :50])
    fig, ax = plt.subplots(figsize=(4, 4))
    leiden = adata.obs["leiden"].astype("category")
    leiden_colors = categorical_colors(len(leiden.cat.categories))
    for index, category in enumerate(leiden.cat.categories):
        selected = np.asarray(leiden == category)
        ax.scatter(
            tsne_embedding[selected, 0], tsne_embedding[selected, 1],
            color=leiden_colors[index], label=str(category), s=2,
            alpha=0.7, rasterized=True,
        )
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_title("t-SNE (perplexity=30, auto learning rate, PCA init)")
    ax.legend(title="Leiden", bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, args.output_dir / "tsne.pdf")

    phate_embedding = phate.PHATE(
        knn=10, decay=40, t="auto", n_jobs=-1, random_state=42
    ).fit_transform(scores[:, :50])
    fig, ax = plt.subplots(figsize=(4, 4))
    if "pseudotime" in adata.obs:
        points = ax.scatter(
            phate_embedding[:, 0], phate_embedding[:, 1],
            c=adata.obs["pseudotime"], cmap="viridis", s=3, alpha=0.7,
            rasterized=True,
        )
        fig.colorbar(points, ax=ax, label="pseudotime")
    else:
        ax.scatter(phate_embedding[:, 0], phate_embedding[:, 1], s=3, alpha=0.7)
    ax.set_xlabel("PHATE 1")
    ax.set_ylabel("PHATE 2")
    ax.set_title("PHATE")
    save_figure(fig, args.output_dir / "phate.pdf")

    retention = neighbor_retention(scores[:, :50], adata.obsm["X_umap"], k=15)
    caption = (
        f"UMAP (n_neighbors=30, min_dist=0.3, random_state=42) of "
        f"N={adata.n_obs} cells. Colored by Leiden cluster (resolution=0.5). "
        f"The 2D map retained {retention:.1%} of 15-nearest-neighbor memberships; "
        "distances between clusters and within-cluster density are not biological signals."
    )
    print(caption)


if __name__ == "__main__":
    main()
