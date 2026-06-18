from markoff_graph import markoff


def test_regression_p31_two_components():
    G = markoff(4, 4, -2, -4, 31)
    assert G.roots() == {(0, 0, 19): 450, (0, 1, 30): 512}


def test_components_matches_roots():
    G = markoff(4, 4, -2, -4, 31)
    assert {r: len(nodes) for r, nodes in G.components().items()} == G.roots()


def test_component_returns_one_root_component():
    G = markoff(4, 4, -2, -4, 31)
    assert len(G.component((0, 0, 19))) == 450
    assert len(G.component((0, 1, 30))) == 512
    assert G.component((0, 0, 19)) == G.components()[(0, 0, 19)]


def test_component_rejects_non_root():
    G = markoff(4, 4, -2, -4, 31)
    try:
        G.component((0, 0, 0))
    except KeyError:
        pass
    else:
        raise AssertionError("expected KeyError")
