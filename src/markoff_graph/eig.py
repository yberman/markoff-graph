"""Optional SciPy eigenvalue helpers."""

from __future__ import annotations

from typing import Tuple

Triple = Tuple[int, int, int]


def _deterministic_v0(triples):
    import numpy as np

    vals = np.empty(len(triples), dtype=np.float64)
    for i, (x, y, z) in enumerate(triples):
        h = 2166136261
        for t in (x, y, z, i):
            h = ((h ^ int(t)) * 16777619) & 0xFFFFFFFF
        vals[i] = float(h & 0xFFFF) - 32768.0
    vals -= vals.mean()
    norm = np.linalg.norm(vals)
    if norm == 0.0:
        vals = np.arange(len(triples), dtype=np.float64)
        vals -= vals.mean()
        norm = np.linalg.norm(vals)
    if norm != 0.0:
        vals /= norm
    return vals


def _normalize_sign(v):
    import numpy as np

    if len(v) == 0:
        return v
    j = int(np.argmax(np.abs(v)))
    if v[j] < 0:
        v = -v
    return v


def component_csr(G, root: Triple):
    """Return ``(A, triples)`` for one component as a SciPy CSR matrix.

    The CSR arrays are built in C as ``double* data``, ``int* indices``, and
    ``int* indptr``.  They are copied into NumPy arrays before the temporary C
    CSR object is freed, then passed directly to ``scipy.sparse.csr_matrix``.
    """
    try:
        from scipy.sparse import csr_matrix
    except ImportError as exc:
        raise ImportError(
            'component_csr requires scipy. Install it with: python -m pip install "markoff-graph[eig]"'
        ) from exc

    data, indices, indptr, _global_indices, triples = G._component_csr_arrays(root)
    n = len(triples)
    return csr_matrix((data, indices, indptr), shape=(n, n)), triples


def _component_csr_full(G, root: Triple):
    try:
        from scipy.sparse import csr_matrix
    except ImportError as exc:
        raise ImportError(
            'eig requires scipy. Install it with: python -m pip install "markoff-graph[eig]"'
        ) from exc

    data, indices, indptr, global_indices, triples = G._component_csr_arrays(root)
    n = len(triples)
    return csr_matrix((data, indices, indptr), shape=(n, n)), global_indices, triples


def eig(G, root: Triple, *, tol: float = 0.0, maxiter=None):
    """Return the second adjacency eigenpair for one connected component.

    The return value is ``(eigenvalue, vector)``, where ``vector`` is a dict
    mapping each triple in the component to its eigenvector coordinate.

    The component adjacency matrix is 3-regular with constant eigenvector of
    eigenvalue 3.  This computes the largest eigenvalue on the subspace
    orthogonal to constants by applying ``P A P`` with SciPy's Lanczos solver.
    """
    try:
        import numpy as np
        from scipy.sparse.linalg import LinearOperator, eigsh
    except ImportError as exc:
        raise ImportError(
            'eig requires numpy and scipy. Install it with: python -m pip install "markoff-graph[eig]"'
        ) from exc

    A, global_indices, triples = _component_csr_full(G, root)
    n = A.shape[0]

    if n == 0:
        return 0.0, {}
    if n == 1:
        vector = np.array([0.0], dtype=np.float64)
        G._store_eig_values(root, 0.0, global_indices, vector)
        return 0.0, {triples[0]: 0.0}

    def project(x):
        return x - np.mean(x)

    def matvec(x):
        y = project(x)
        z = A @ y
        return project(z)

    op = LinearOperator((n, n), matvec=matvec, dtype=np.float64)
    v0 = _deterministic_v0(triples)
    values, vectors = eigsh(op, k=1, which="LA", v0=v0, tol=tol, maxiter=maxiter)

    eigenvalue = float(values[0])
    vector = np.asarray(vectors[:, 0], dtype=np.float64)
    vector = project(vector)
    norm = np.linalg.norm(vector)
    if norm != 0.0:
        vector = vector / norm
    vector = _normalize_sign(vector)

    G._store_eig_values(root, eigenvalue, global_indices, vector)
    return eigenvalue, {triple: float(value) for triple, value in zip(triples, vector)}
