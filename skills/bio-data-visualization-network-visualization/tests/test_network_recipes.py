"""Focused regressions for the audited network-visualization failures.

Usage: python tests/test_network_recipes.py
"""

from __future__ import annotations

import importlib.util
import tempfile
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

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

    large = nx.path_graph([f"G{i:03d}" for i in range(100)])
    labels = module.adaptive_labels(large, genes_of_interest={"G099"})
    assert len(labels) == 16
    assert labels["G099"] == "G099"
    calls = []

    def capture_labels(*_args, **kwargs):
        calls.append(kwargs.get("labels"))

    with tempfile.TemporaryDirectory() as directory:
        with patch.object(module.nx, "draw_networkx_labels", capture_labels):
            module.render_networks(large, Path(directory), genes_of_interest={"G099"})
    assert len(calls) == 4
    assert all(call == labels for call in calls)


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
    assert html.count("<h1>") == 1
    assert html.count("Biological Network") == 1


def test_signed_edges_and_cytoscape_labels() -> None:
    directed = load_module("directed_network", "scripts/directed_network.py")
    graph = nx.DiGraph()
    graph.add_edge("TF", "target", sign="activation")
    try:
        directed.validate_sign_values(graph, "sign")
    except ValueError as error:
        assert "activation" in str(error)
    else:
        raise AssertionError("unknown regulatory sign was accepted")

    cytoscape = load_module("cytoscape_automation", "examples/cytoscape_automation.py")
    ppi = nx.path_graph([f"G{i:03d}" for i in range(100)])
    cytoscape.prepare_style_attributes(ppi)
    labels = nx.get_node_attributes(ppi, "display_label")
    assert len(labels) == 100
    assert sum(bool(label) for label in labels.values()) == 15

    calls = []
    with (
        patch.object(cytoscape.p4c, "create_network_from_networkx", return_value=7),
        patch.object(cytoscape.p4c, "layout_network", side_effect=lambda *_: calls.append("layout")),
        patch.object(cytoscape.p4c, "fit_content", side_effect=lambda *_args, **_kwargs: calls.append("fit")),
    ):
        cytoscape.send_network_to_cytoscape(ppi)
    assert calls == ["layout", "fit"]

    with (
        patch.object(cytoscape.p4c, "create_visual_style"),
        patch.object(cytoscape.p4c, "set_node_size_mapping"),
        patch.object(cytoscape.p4c, "set_node_color_mapping"),
        patch.object(cytoscape.p4c, "set_edge_line_width_mapping"),
        patch.object(cytoscape.p4c, "set_node_shape_mapping"),
        patch.object(cytoscape.p4c, "set_node_label_mapping") as label_mapping,
        patch.object(cytoscape.p4c, "set_visual_style"),
    ):
        cytoscape.apply_ppi_style()
    label_mapping.assert_called_once_with("display_label", style_name="PPIStyle")


def main() -> None:
    test_static_helpers()
    test_pyvis_preserves_attributes_and_direction()
    test_signed_edges_and_cytoscape_labels()
    print("PASS: adaptive labels, sign validation, PyVis heading, Cytoscape labels")


if __name__ == "__main__":
    main()
