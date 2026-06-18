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


def test_capacity_error_is_reported():
    try:
        markoff(4, 4, -2, -4, 31, capacity=1)
    except RuntimeError as exc:
        assert "too many solutions" in str(exc)
    else:
        raise AssertionError("expected RuntimeError")


def test_prime_bound_is_not_old_1621_bound():
    try:
        markoff(4, 4, -2, -4, 1627, capacity=1)
    except RuntimeError as exc:
        assert "too many solutions" in str(exc)
    else:
        raise AssertionError("expected RuntimeError")


def test_rejects_prime_above_uint32_node_bound():
    try:
        markoff(4, 4, -2, -4, 46349)
    except ValueError as exc:
        assert "2*p*p must fit in uint32" in str(exc)
    else:
        raise AssertionError("expected ValueError")
