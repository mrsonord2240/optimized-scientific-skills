"""Interactive HTML network visualization with PyVis.

Inputs: optional GraphML; otherwise a reproducible PPI-like demo graph.
Usage: ``python examples/interactive_network.py --graphml network.graphml --output-dir out/interactive``
"""

from __future__ import annotations

import argparse
import html
import json
import re
from pathlib import Path

import networkx as nx
import numpy as np
from networkx.algorithms.community import greedy_modularity_communities
from pyvis.network import Network


PALETTE = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
]


def save_with_single_heading(net: Network, output: Path, title: str) -> Path:
    """Save PyVis HTML and replace template-duplicated headings with one escaped title."""
    output.parent.mkdir(parents=True, exist_ok=True)
    net.save_graph(str(output))
    document = output.read_text(encoding="utf-8")
    document = re.sub(
        r"\s*<center>\s*<h1>.*?</h1>\s*</center>\s*",
        "\n",
        document,
        flags=re.DOTALL,
    )
    heading = f"<center><h1>{html.escape(title)}</h1></center>"
    if "<body>" not in document:
        raise RuntimeError("PyVis template did not emit a body element")
    document = document.replace("<body>", f"<body>\n{heading}", 1)
    output.write_text(document, encoding="utf-8", newline="\n")
    return output


def build_demo_graph() -> nx.Graph:
    """Build a deterministic PPI-like graph with confidence weights."""
    rng = np.random.default_rng(42)
    genes = [f"Gene{i}" for i in range(1, 41)]
    graph = nx.barabasi_albert_graph(40, 2, seed=42)
    graph = nx.relabel_nodes(graph, {i: genes[i] for i in range(40)})
    for u, v in graph.edges():
        graph[u][v]["weight"] = float(rng.uniform(0.4, 1.0))
        graph[u][v]["score"] = int(graph[u][v]["weight"] * 1000)
    return graph


def create_interactive_network(
    graph: nx.Graph, output: Path, title: str = "Biological Network"
) -> Path:
    """Write a basic network without allowing PyVis to mutate the caller's graph."""
    net = Network(
        height="700px", width="100%", bgcolor="white", font_color="black",
        heading="", directed=graph.is_directed(), cdn_resources="in_line",
    )
    net.from_nx(graph.copy())
    net.toggle_physics(True)
    net.show_buttons(filter_=["physics"])
    return save_with_single_heading(net, output, title)


def create_styled_network(
    graph: nx.Graph, output: Path, title: str = "Styled Network"
) -> Path:
    """Write an interactive network with degree sizing and community coloring."""
    community_graph = graph.to_undirected() if graph.is_directed() else graph
    communities = list(greedy_modularity_communities(community_graph))
    membership = {
        node: index for index, community in enumerate(communities) for node in community
    }
    degrees = dict(graph.degree())
    net = Network(
        height="700px", width="100%", bgcolor="white", font_color="black",
        heading="", directed=graph.is_directed(), cdn_resources="in_line",
    )
    for node in graph.nodes():
        community = membership.get(node, 0)
        tooltip = f"{node}\nDegree: {degrees[node]}\nCommunity: {community + 1}"
        net.add_node(
            node, size=10 + degrees[node] * 5,
            color=PALETTE[community % len(PALETTE)], title=tooltip, label=node,
        )
    for u, v, data in graph.edges(data=True):
        weight = float(data.get("weight", 0.5))
        score = data.get("score", 500)
        net.add_edge(u, v, width=weight * 3, title=f"{u} -- {v}\nScore: {score}")
    net.set_options(json.dumps({
        "physics": {
            "forceAtlas2Based": {
                "gravitationalConstant": -50, "centralGravity": 0.01,
                "springLength": 100, "springConstant": 0.08,
            },
            "solver": "forceAtlas2Based", "stabilization": {"iterations": 150},
        },
        "interaction": {"hover": True, "navigationButtons": True, "keyboard": True},
    }))
    return save_with_single_heading(net, output, title)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graphml", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("out/interactive"))
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    graph = nx.read_graphml(args.graphml) if args.graphml else build_demo_graph()
    basic = create_interactive_network(
        graph, output_dir / "network_basic.html", "Basic Interactive Network"
    )
    styled = create_styled_network(
        graph, output_dir / "network_styled.html", "PPI Network with Communities"
    )
    communities = list(greedy_modularity_communities(graph))
    print(f"Basic network saved: {basic}")
    print(f"Styled network saved: {styled}")
    print(f"Network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    print(f"Communities detected: {len(communities)}")


if __name__ == "__main__":
    main()
