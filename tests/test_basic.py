import math

import pytest

from markoff_graph import MarkoffGraph


def test_regression_p31_two_components():
    G = MarkoffGraph(4, 4, -2, -4, 31)
    assert G.roots() == {(0, 0, 19): 450, (0, 1, 30): 512}


def test_regression_p43_two_components():
    G = MarkoffGraph(4, 4, -2, -4, 43)
    assert G.roots() == {(0, 0, 29): 882, (0, 2, 0): 968}


def test_nodes_count_matches_roots():
    G = MarkoffGraph(4, 4, -2, -4, 31)
    assert len(list(G.nodes())) == sum(G.roots().values())


def test_component_returns_one_root_component():
    G = MarkoffGraph(4, 4, -2, -4, 31)
    assert len(G.component((0, 0, 19))) == 450
    assert len(G.component((0, 1, 30))) == 512


def test_component_rejects_non_root():
    G = MarkoffGraph(4, 4, -2, -4, 31)
    with pytest.raises(KeyError):
        G.component((0, 0, 0))


def test_public_methods_are_small():
    G = MarkoffGraph(4, 4, -2, -4, 31)
    assert hasattr(G, "nodes")
    assert hasattr(G, "roots")
    assert hasattr(G, "component")
    assert hasattr(G, "eig")
    assert not hasattr(G, "edges")
    assert not hasattr(G, "neighbors")
    assert not hasattr(G, "components")


def test_rejects_prime_above_uint32_node_bound():
    with pytest.raises(ValueError, match="2\\*p\\*p must fit in uint32"):
        MarkoffGraph(4, 4, -2, -4, 46349)


def test_context_manager_closes_graph():
    with MarkoffGraph(4, 4, -2, -4, 31) as G:
        assert len(list(G.nodes())) == 962
    with pytest.raises(RuntimeError, match="closed"):
        list(G.nodes())


def test_nx_export_if_networkx_is_available():
    pytest.importorskip("networkx")
    from markoff_graph.nx_export import to_multigraph

    G = MarkoffGraph(4, 4, -2, -4, 31)
    H = to_multigraph(G)
    assert H.number_of_nodes() == 962
    assert H.number_of_nodes() == sum(G.roots().values())

    H0 = to_multigraph(G, root=(0, 0, 19))
    assert H0.number_of_nodes() == 450


def test_native_component_data_matches_roots():
    G = MarkoffGraph(4, 4, -2, -4, 31)
    seen = {}
    for i in range(G._component_count):
        c = G._component_at(i)
        root_node = G._node_at(int(c.root_index))
        root = (int(root_node.x), int(root_node.y), int(root_node.z))
        assert c.root.contents.x == root_node.x
        assert c.root.contents.y == root_node.y
        assert c.root.contents.z == root_node.z
        assert int(c.size) == G.roots()[root]
        seen[root] = int(c.size)
    assert seen == G.roots()


def test_csr_export_if_scipy_is_available():
    pytest.importorskip("scipy")
    from markoff_graph.eig import component_csr

    G = MarkoffGraph(4, 4, -2, -4, 31)
    A, triples = component_csr(G, (0, 0, 19))
    assert A.shape == (450, 450)
    assert len(triples) == 450
    assert A.nnz == 3 * 450
    row_sums = A.sum(axis=1).A.ravel()
    assert row_sums.min() == 3.0
    assert row_sums.max() == 3.0


def test_scipy_eig_p43_if_available():
    pytest.importorskip("scipy")

    G = MarkoffGraph(4, 4, -2, -4, 43)
    value, vector = G.eig((0, 0, 29))
    assert math.isclose(value, 2.92233598310223, rel_tol=0.0, abs_tol=1e-10)
    assert set(vector) == G.component((0, 0, 29))
    assert abs(sum(vector.values())) < 1e-10
    assert abs(sum(x * x for x in vector.values()) - 1.0) < 1e-10

    root_index = G._component_root_index((0, 0, 29))
    stored = [
        float(G._node_at(j).eigenvector)
        for j in range(G._node_count)
        if int(G._node_at(j).root) == root_index
    ]
    assert abs(sum(stored)) < 1e-10
    assert abs(sum(x * x for x in stored) - 1.0) < 1e-10


def test_scipy_eig_p41_if_available():
    pytest.importorskip("scipy")

    G = MarkoffGraph(2, 1, 0, -5, 41)
    assert G.roots() == {(0, 0, 6): 1723}

    value, vector = G.eig((0, 0, 6))
    assert math.isclose(value, 2.876744056040619, rel_tol=0.0, abs_tol=1e-10)
    assert set(vector) == G.component((0, 0, 6))
    assert abs(sum(vector.values())) < 1e-10
    assert abs(sum(x * x for x in vector.values()) - 1.0) < 1e-10
