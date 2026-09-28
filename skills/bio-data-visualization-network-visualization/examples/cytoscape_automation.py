"""Automate a publication-quality PPI figure in a running Cytoscape session.

Inputs: optional GraphML; otherwise a reproducible PPI-like demo graph.
Usage: ``python examples/cytoscape_automation.py --graphml network.graphml --output-dir out/cytoscape``
"""

from __future__ import annotations

import argparse
import traceback
from pathlib import Path

import networkx as nx
import numpy as np
import py4cytoscape as p4c


def build_demo_graph() -> nx.Graph:
    """Build a deterministic PPI-like graph with Cytoscape style attributes."""
    rng = np.random.default_rng(42)
    genes = [
        "TP53", "MDM2", "BRCA1", "ATM", "CHEK2", "CDK2", "RB1", "CDKN1A",
        "BAX", "BCL2", "CASP3", "CASP9", "BID", "AKT1", "MAPK1",
    ]
    graph = nx.barabasi_albert_graph(len(genes), 2, seed=42)
    graph = nx.relabel_nodes(graph, {i: genes[i] for i in range(len(genes))})
    kinases = {"ATM", "CHEK2", "CDK2", "AKT1", "MAPK1"}
    nx.set_node_attributes(graph, dict(graph.degree()), "degree")
    nx.set_node_attributes(
        graph,
        {node: "kinase" if node in kinases else "other" for node in graph},
        "gene_type",
    )
    for u, v in graph.edges():
        graph[u][v]["score"] = int(rng.integers(400, 1000))
        graph[u][v]["weight"] = graph[u][v]["score"] / 1000
    return graph


def prepare_style_attributes(graph: nx.Graph) -> nx.Graph:
    """Populate every column referenced by the Cytoscape visual mappings."""
    nx.set_node_attributes(graph, dict(graph.degree()), "degree")
    for node in graph:
        graph.nodes[node].setdefault("gene_type", "other")
    for _, _, data in graph.edges(data=True):
        if "score" not in data:
            data["score"] = int(float(data.get("weight", 1.0)) * 1000)
        data.setdefault("weight", float(data["score"]) / 1000)
    return graph


def send_network_to_cytoscape(graph: nx.Graph, title: str = "PPI Network") -> int:
    """Send a NetworkX graph to Cytoscape and apply a force-directed layout."""
    suid = p4c.create_network_from_networkx(graph, title=title)
    p4c.layout_network("force-directed")
    print(f"Network created in Cytoscape (SUID: {suid})")
    return suid


def apply_ppi_style() -> str:
    """Create and apply a PPI-appropriate visual style."""
    style_name = "PPIStyle"
    defaults = {
        "NODE_SHAPE": "ELLIPSE",
        "NODE_FILL_COLOR": "#4DBBD5",
        "NODE_BORDER_WIDTH": 1.0,
        "NODE_BORDER_PAINT": "#000000",
        "NODE_LABEL_FONT_SIZE": 10,
        "EDGE_STROKE_UNSELECTED_PAINT": "#999999",
        "EDGE_TARGET_ARROW_SHAPE": "NONE",
        "NETWORK_BACKGROUND_PAINT": "#FFFFFF",
    }
    p4c.create_visual_style(style_name, defaults=defaults)
    p4c.set_node_size_mapping(
        "degree", [1, 5, 15], [30, 60, 120], mapping_type="c", style_name=style_name
    )
    p4c.set_node_color_mapping(
        "degree", [1, 5, 15], ["#FFFFCC", "#FD8D3C", "#BD0026"],
        mapping_type="c", style_name=style_name,
    )
    p4c.set_edge_line_width_mapping(
        "score", [400, 700, 1000], [1, 2, 4],
        mapping_type="c", style_name=style_name,
    )
    # Shape mapping is always discrete in py4cytoscape 1.13 and has no mapping_type.
    p4c.set_node_shape_mapping(
        "gene_type", ["kinase", "other"], ["DIAMOND", "ELLIPSE"],
        style_name=style_name,
    )
    p4c.set_visual_style(style_name)
    print(f"Applied style: {style_name}")
    return style_name


def export_figure(filename: Path, resolution: int = 300) -> Path:
    """Export the current network view to an absolute, overwriteable path."""
    path = filename.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        p4c.export_image(str(path), type="PDF", overwrite_file=True)
    elif suffix == ".svg":
        p4c.export_image(str(path), type="SVG", overwrite_file=True)
    else:
        p4c.export_image(
            str(path), type="PNG", resolution=resolution, overwrite_file=True
        )
    print(f"Exported: {path}")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graphml", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("out/cytoscape"))
    args = parser.parse_args()
    p4c.cytoscape_ping()
    print("Cytoscape is running.")
    graph = nx.read_graphml(args.graphml) if args.graphml else build_demo_graph()
    send_network_to_cytoscape(prepare_style_attributes(graph))
    apply_ppi_style()
    export_figure(args.output_dir / "ppi_network.pdf")
    export_figure(args.output_dir / "ppi_network.png", resolution=300)
    print("Network visualization complete.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
