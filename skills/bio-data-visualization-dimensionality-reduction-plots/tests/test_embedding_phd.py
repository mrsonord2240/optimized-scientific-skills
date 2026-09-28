"""Focused regressions for examples/embedding_phd.py."""

from __future__ import annotations

import importlib.util
from itertools import combinations
from pathlib import Path
import unittest

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EXAMPLE = Path(__file__).parents[1] / "examples" / "embedding_phd.py"
SPEC = importlib.util.spec_from_file_location("embedding_phd", EXAMPLE)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class EmbeddingExampleTests(unittest.TestCase):
    def test_example_selects_hvgs_before_normalization(self) -> None:
        source = EXAMPLE.read_text(encoding="utf-8")
        self.assertLess(
            source.index("sc.pp.highly_variable_genes"),
            source.index("sc.pp.normalize_total"),
        )
        for filename in (
            "pca.pdf", "scree.pdf", "umap_clusters.pdf", "tsne.pdf", "phate.pdf"
        ):
            self.assertIn(filename, source)

    def test_neighbor_retention_identity_is_one(self) -> None:
        rng = np.random.default_rng(42)
        values = rng.normal(size=(30, 2))
        self.assertEqual(MODULE.neighbor_retention(values, values.copy(), k=5), 1.0)

    def test_twelve_categories_receive_unique_rgba_values(self) -> None:
        colors = MODULE.categorical_colors(12)
        self.assertEqual(colors.shape, (12, 4))
        self.assertEqual(len({tuple(color) for color in colors}), 12)
        with self.assertRaisesRegex(ValueError, "facet"):
            MODULE.categorical_colors(21)

    def test_loading_labels_do_not_overlap(self) -> None:
        fig, ax = plt.subplots(figsize=(5, 4))
        endpoints = np.zeros((5, 2))
        labels = [f"LONG_GENE_LABEL_{index}" for index in range(5)]
        annotations = MODULE.annotate_loading_labels(ax, endpoints, labels)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [annotation.get_window_extent(renderer) for annotation in annotations]
        self.assertTrue(
            all(not left.overlaps(right) for left, right in combinations(boxes, 2))
        )
        self.assertEqual(
            [annotation.get_position() for annotation in annotations],
            [(6, 0), (6, -12), (6, 12), (6, -24), (6, 24)],
        )
        plt.close(fig)

    def test_validate_input_accepts_raw_counts(self) -> None:
        rng = np.random.default_rng(42)
        adata = ad.AnnData(
            rng.poisson(1, size=(100, 2_000)),
            obs=pd.DataFrame({"condition": ["control"] * 100}),
        )
        MODULE.validate_input(adata)

    def test_validate_input_rejects_log_values(self) -> None:
        adata = ad.AnnData(
            np.full((100, 2_000), 0.5),
            obs=pd.DataFrame({"condition": ["control"] * 100}),
        )
        with self.assertRaisesRegex(ValueError, "raw"):
            MODULE.validate_input(adata)

    def test_skill_references_resolve_and_main_is_short(self) -> None:
        skill_dir = EXAMPLE.parents[1]
        skill = (skill_dir / "SKILL.md").read_text(encoding="utf-8")
        self.assertLessEqual(len(skill.splitlines()), 300)
        self.assertIn(
            "name: bio-data-visualization-dimensionality-reduction-plots", skill
        )
        for relative in (
            "references/method-recipes.md",
            "references/neighborhood-validation.md",
            "examples/embedding_phd.py",
        ):
            self.assertTrue((skill_dir / relative).is_file(), relative)
        self.assertEqual(skill.count("```"), 0)


if __name__ == "__main__":
    unittest.main()
