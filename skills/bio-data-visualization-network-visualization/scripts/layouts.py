"""Render reproducible NetworkX layouts from a GraphML network.

Inputs: GraphML; output directory; optional comma-separated layout names.
Usage: python scripts/layouts.py network.graphml --output-dir out/layouts
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from networkx.algorithms.community import greedy_modularity_communities


def compute_layout(graph: nx.Graph, name: str) -> dict:
    """Compute one supported layout with deterministic stochastic initialization."""
    if name == "spring":
        return nx.spring_layout(graph, seed=42)
    if name == "kamada-kawai":
        return nx.kamada_kawai_layout(graph)
    if name == "circular":
        return nx.circular_layout(graph)
    if name == "spectral":
        return nx.spectral_layout(graph)
    if name == "forceatlas2":
        return nx.forceatlas2_layout(graph, seed=42)
    raise ValueError(f"Unsupported layout: {name}")


def render_layout(graph: nx.Graph, name: str, output: Path) -> Path:
    """Render one layout with degree sizing, community color, and adaptive labels."""
    positions = compute_layout(graph, name)
    if set(positions) != set(graph):
        raise RuntimeError(f"{name} omitted graph nodes")
    xy = np.asarray(list(positions.values()), dtype=float)
    if not np.isfinite(xy).all():
        raise RuntimeError(f"{name} produced non-finite coordinates")

    community_graph = graph.to_undirected() if graph.is_directed() else graph
    communities = list(greedy_modularity_communities(community_graph))
    membership = {
        node: index for index, community in enumerate(communities) for node in community
    }
    degrees = dict(graph.degree())
    label_count = min(len(graph), 30, max(15, math.ceil(0.02 * len(graph))))
    labelled = sorted(degrees, key=degrees.get, reverse=True)[:label_count]

    fig, ax = plt.subplots(figsize=(10, 8))
    nx.draw_networkx_edges(graph, positions, alpha=0.25, width=0.7, ax=ax)
    nx.draw_networkx_nodes(
        graph,
        positions,
        node_size=[80 + 30 * degrees[node] for node in graph],
        node_color=[membership[node] for node in graph],
        cmap="tab20",
        edgecolors="black",
        linewidths=0.4,
        ax=ax,
    )
    nx.draw_networkx_labels(
        graph,
        positions,
        labels={node: str(node) for node in labelled},
        font_size=7,
        ax=ax,
    )
    ax.set_title(f"{name} layout")
    ax.axis("off")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"{name}: {len(positions)} positions; {output.stat().st_size} bytes")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graphml", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--layouts",
        default="spring,kamada-kawai,circular,spectral,forceatlas2",
        help="Comma-separated: spring, kamada-kawai, circular, spectral, forceatlas2",
    )
    args = parser.parse_args()
    graph = nx.read_graphml(args.graphml)
    if not graph:
        raise ValueError("The input graph is empty")
    for name in args.layouts.split(","):
        render_layout(graph, name.strip(), args.output_dir.resolve() / f"{name.strip()}.png")


if __name__ == "__main__":
    main()
