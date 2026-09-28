"""Focused regressions for the audited network-visualization failures.

Usage: python tests/test_network_recipes.py
"""

from __future__ import annotations

import importlib.util
import tempfile
from copy import deepcopy
from pathlib import Path

import networkx as nx


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {relative_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_static_helpers() -> None:
    module = load_module("network_plots", "examples/network_plots.py")
    unweighted = nx.path_graph(4)
    widths = module.visible_edge_widths(unweighted)
    assert len(widths) == 3
    assert all(0.5 <= width <= 4.0 for width in widths)

    graph = nx.karate_club_graph()
    communities, membership, colors = module.community_style(graph)
    palette = module.plt.colormaps["Set2"].resampled(len(communities))
    expected = [palette(membership[node]) for node in graph]
    assert colors == expected


def test_pyvis_preserves_attributes_and_direction() -> None:
    module = load_module("interactive_network", "examples/interactive_network.py")
    graph = nx.DiGraph()
    graph.add_edge("TF", "target", weight=0.75, score=750)
    before = deepcopy(dict(graph.edges))
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "directed.html"
        module.create_interactive_network(graph, output)
        html = output.read_text(encoding="utf-8")
    assert dict(graph.edges) == before
    assert graph["TF"]["target"]["weight"] == 0.75
    assert '"arrows": "to"' in html or '"arrows": "to"'.replace(" ", "") in html.replace(" ", "")


def main() -> None:
    test_static_helpers()
    test_pyvis_preserves_attributes_and_direction()
    print("PASS: static mapping, unweighted widths, PyVis immutability, directed arrows")


if __name__ == "__main__":
    main()
