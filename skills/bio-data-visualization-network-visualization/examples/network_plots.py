"""Static network visualization with NetworkX and matplotlib.

Inputs: optional GraphML; otherwise a reproducible PPI-like demo graph.
Usage: ``python examples/network_plots.py --graphml network.graphml --output-dir out/static``
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from networkx.algorithms.community import greedy_modularity_communities


def build_demo_graph() -> nx.Graph:
    """Build a deterministic PPI-like graph with STRING-style scores."""
    rng = np.random.default_rng(42)
    genes = [f"Gene{i}" for i in range(1, 31)]
    graph = nx.barabasi_albert_graph(30, 2, seed=42)
    graph = nx.relabel_nodes(graph, {i: genes[i] for i in range(30)})
    for u, v in graph.edges():
        graph[u][v]["score"] = int(rng.integers(400, 1000))
        graph[u][v]["weight"] = graph[u][v]["score"] / 1000
    return graph


def visible_edge_widths(
    graph: nx.Graph, minimum: float = 0.5, maximum: float = 4.0
) -> list[float]:
    """Scale optional edge weights into a visible width range."""
    raw = np.asarray(
        [float(data.get("weight", 1.0)) for _, _, data in graph.edges(data=True)],
        dtype=float,
    )
    if raw.size == 0:
        return []
    span = float(np.ptp(raw))
    if span == 0:
        return [float((minimum + maximum) / 2)] * len(raw)
    return (minimum + (raw - raw.min()) * (maximum - minimum) / span).tolist()


def community_style(
    graph: nx.Graph,
) -> tuple[list[set[str]], dict[str, int], list[tuple[float, float, float, float]]]:
    """Return communities, memberships, and node colors from one shared mapping."""
    communities = [set(group) for group in greedy_modularity_communities(graph)]
    membership = {
        node: index for index, community in enumerate(communities) for node in community
    }
    palette = plt.colormaps["Set2"].resampled(max(1, len(communities)))
    colors = [palette(membership[node]) for node in graph.nodes()]
    return communities, membership, colors


def adaptive_labels(
    graph: nx.Graph, genes_of_interest: Iterable[str] = ()
) -> dict[str, str]:
    """Return capped degree-ranked labels plus requested nodes present in the graph."""
    degrees = dict(graph.degree())
    label_count = min(len(graph), 30, max(15, math.ceil(0.02 * len(graph))))
    ranked = sorted(graph, key=lambda node: (-degrees[node], str(node)))[:label_count]
    requested = [node for node in genes_of_interest if node in graph]
    selected = dict.fromkeys([*ranked, *requested])
    return {node: str(node) for node in selected}


def render_networks(
    graph: nx.Graph, output_dir: Path, genes_of_interest: Iterable[str] = ()
) -> list[Path]:
    """Render four static network views and return their paths."""
    output_dir.mkdir(parents=True, exist_ok=True)
    pos = nx.spring_layout(graph, seed=42, k=1.5)
    degrees = dict(graph.degree())
    labels = adaptive_labels(graph, genes_of_interest)
    node_sizes = [100 + degrees[node] * 150 for node in graph.nodes()]
    edge_widths = visible_edge_widths(graph)
    outputs: list[Path] = []

    fig, ax = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_edges(
        graph, pos, width=edge_widths, alpha=0.3, edge_color="gray", ax=ax
    )
    nodes = nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=node_sizes,
        node_color=[degrees[node] for node in graph.nodes()],
        cmap=plt.cm.YlOrRd,
        edgecolors="black",
        linewidths=0.5,
        ax=ax,
    )
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=7, ax=ax)
    fig.colorbar(nodes, ax=ax, label="Degree", shrink=0.8)
    ax.set_title(
        f"PPI Network ({graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges)"
    )
    ax.axis("off")
    outputs.append(output_dir / "network_degree.png")
    fig.savefig(outputs[-1], dpi=300, bbox_inches="tight")
    plt.close(fig)

    communities, _membership, community_colors = community_style(graph)
    palette = plt.colormaps["Set2"].resampled(max(1, len(communities)))
    fig, ax = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_edges(graph, pos, alpha=0.2, edge_color="gray", ax=ax)
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=node_sizes,
        node_color=community_colors,
        edgecolors="black",
        linewidths=0.5,
        ax=ax,
    )
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=7, ax=ax)
    for index, community in enumerate(communities):
        ax.scatter(
            [], [], c=[palette(index)], s=80,
            label=f"Community {index + 1} ({len(community)} nodes)",
        )
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.set_title(f"Community Structure ({len(communities)} communities)")
    ax.axis("off")
    outputs.append(output_dir / "network_communities.png")
    fig.savefig(outputs[-1], dpi=300, bbox_inches="tight")
    plt.close(fig)

    hub_genes = set(sorted(degrees, key=degrees.get, reverse=True)[:5])
    fig, ax = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_edges(graph, pos, alpha=0.15, edge_color="gray", ax=ax)
    nx.draw_networkx_nodes(
        graph, pos,
        node_size=[800 if node in hub_genes else 200 for node in graph.nodes()],
        node_color=["#E64B35" if node in hub_genes else "#cccccc" for node in graph.nodes()],
        edgecolors="black", linewidths=0.5, ax=ax,
    )
    nx.draw_networkx_labels(
        graph, pos, labels=labels,
        font_size=10, font_weight="bold", ax=ax,
    )
    ax.set_title("Hub Genes Highlighted (top 5 by degree)")
    ax.axis("off")
    outputs.append(output_dir / "network_hubs.png")
    fig.savefig(outputs[-1], dpi=300, bbox_inches="tight")
    plt.close(fig)

    edge_colors = [
        "#E64B35" if data.get("score", 0) >= 900
        else "#F39B7F" if data.get("score", 0) >= 700
        else "#cccccc"
        for _, _, data in graph.edges(data=True)
    ]
    fig, ax = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_edges(
        graph, pos, edge_color=edge_colors, width=edge_widths, alpha=0.7, ax=ax
    )
    nx.draw_networkx_nodes(
        graph, pos, node_size=node_sizes, node_color="#4DBBD5",
        edgecolors="black", linewidths=0.5, ax=ax,
    )
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=7, ax=ax)
    ax.scatter([], [], c="#E64B35", s=60, label="Score >= 900")
    ax.scatter([], [], c="#F39B7F", s=60, label="Score >= 700")
    ax.scatter([], [], c="#cccccc", s=60, label="Score < 700")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9, title="Confidence")
    ax.set_title("PPI Network Colored by Interaction Confidence")
    ax.axis("off")
    outputs.append(output_dir / "network_confidence.png")
    fig.savefig(outputs[-1], dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    print(f"Density: {nx.density(graph):.3f}")
    print(f"Communities: {len(communities)}")
    print(f"Avg clustering coefficient: {nx.average_clustering(graph):.3f}")
    degree_df = pd.DataFrame(
        sorted(graph.degree(), key=lambda item: item[1], reverse=True),
        columns=["gene", "degree"],
    )
    print("\nTop 10 by degree:")
    print(degree_df.head(10).to_string(index=False))
    print("\nPlots saved:", ", ".join(str(path) for path in outputs))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graphml", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("out/static"))
    parser.add_argument(
        "--label", action="append", default=[], metavar="NODE",
        help="also label this node (repeat for multiple genes of interest)",
    )
    args = parser.parse_args()
    graph = nx.read_graphml(args.graphml) if args.graphml else build_demo_graph()
    render_networks(graph, args.output_dir.resolve(), args.label)


if __name__ == "__main__":
    main()
