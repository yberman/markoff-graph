"""Optional NetworkX export helpers."""

from __future__ import annotations

from typing import Optional, Tuple

Triple = Tuple[int, int, int]


def to_multigraph(G, root: Optional[Triple] = None):
    """Convert a MarkoffGraph, or one component, to a NetworkX MultiGraph.

    NetworkX is an optional dependency. Install it with:

        python -m pip install "markoff-graph[nx]"
    """
    try:
        import networkx as nx
    except ImportError as exc:
        raise ImportError(
            "to_multigraph requires networkx. Install it with: "
            'python -m pip install "markoff-graph[nx]"'
        ) from exc

    if root is None:
        nodes = set(G.nodes())
    else:
        nodes = set(G.component(root))

    H = nx.MultiGraph()
    H.add_nodes_from(nodes)

    A = G._A
    B = G._B
    C = G._C
    p = G._prime

    for x, y, z in nodes:
        u = (x, y, z)
        neighbors = (
            ("sigma_x", ((y * z + A - x) % p, y, z)),
            ("sigma_y", (x, (x * z + B - y) % p, z)),
            ("sigma_z", (x, y, (x * y + C - z) % p)),
        )
        for move, v in neighbors:
            if v in nodes and u <= v:
                H.add_edge(u, v, key=move, move=move)

    return H
