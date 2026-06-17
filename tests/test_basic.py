from markoff_graph import markoff


def test_regression_p31_two_components():
    G = markoff(4, 4, -2, -4, 31)
    assert G.roots() == {(0, 0, 19): 450, (0, 1, 30): 512}


def test_component_matches_roots():
    G = markoff(4, 4, -2, -4, 31)
    assert {r: len(nodes) for r, nodes in G.component().items()} == G.roots()
