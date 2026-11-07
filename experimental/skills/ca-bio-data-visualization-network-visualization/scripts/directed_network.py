"""Render a directed signed regulatory network with Graphviz hierarchy.

Inputs: directed GraphML; edge attribute containing '+' and '-' signs.
Usage: python scripts/directed_network.py grn.graphml --output out/grn.png
Requires: pydot and the Graphviz ``dot`` executable.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D


SIGN_COLORS = {"+": "#D55E00", "-": "#0072B2"}


def validate_sign_values(graph: nx.DiGraph, sign_attribute: str) -> None:
    """Reject explicit sign values that the renderer would otherwise omit."""
    values = {
        str(data.get(sign_attribute, "+")) for _, _, data in graph.edges(data=True)
    }
    unknown = sorted(values - set(SIGN_COLORS))
    if unknown:
        raise ValueError(
            f"Unknown {sign_attribute!r} value(s): {unknown}; "
            "normalize activation/repression to '+'/'-' before rendering"
        )


def render_directed(graph: nx.DiGraph, output: Path, sign_attribute: str) -> Path:
    """Render direction with arrowheads and regulation sign with edge color."""
    if not graph.is_directed():
        raise ValueError("Expected a directed graph")
    validate_sign_values(graph, sign_attribute)
    positions = nx.nx_pydot.graphviz_layout(graph, prog="dot")
    fig, ax = plt.subplots(figsize=(11, 8))
    nx.draw_networkx_nodes(
        graph,
        positions,
        node_size=[140 + 35 * graph.out_degree(node) for node in graph],
        node_color="#E5E5E5",
        edgecolors="black",
        linewidths=0.6,
        ax=ax,
    )
    for sign, color in SIGN_COLORS.items():
        edges = [
            (u, v)
            for u, v, data in graph.edges(data=True)
            if str(data.get(sign_attribute, "+")) == sign
        ]
        nx.draw_networkx_edges(
            graph,
            positions,
            edgelist=edges,
            edge_color=color,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=12,
            width=1.0,
            alpha=0.75,
            ax=ax,
        )
    regulators = [
        node
        for node in sorted(graph, key=graph.out_degree, reverse=True)
        if graph.out_degree(node) > 0
    ][:15]
    nx.draw_networkx_labels(
        graph,
        positions,
        labels={node: str(node) for node in regulators},
        font_size=7,
        ax=ax,
    )
    ax.legend(
        handles=[
            Line2D([0], [0], color=SIGN_COLORS["+"], label="activation (+)"),
            Line2D([0], [0], color=SIGN_COLORS["-"], label="repression (-)"),
        ],
        loc="upper left",
    )
    ax.axis("off")
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    if output.stat().st_size == 0:
        raise RuntimeError("Renderer created an empty file")
    print(
        f"directed: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges; "
        f"{output.stat().st_size} bytes"
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graphml", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sign-attribute", default="sign")
    args = parser.parse_args()
    graph = nx.read_graphml(args.graphml)
    render_directed(graph, args.output, args.sign_attribute)


if __name__ == "__main__":
    main()
