"""Finite-field Markoff graph components.

Core usage:

    from markoff_graph import MarkoffGraph

    G = MarkoffGraph(4, 4, -2, -4, 31)
    print(G.roots())

The core graph object exposes only:

    G.nodes()          iterator of solution triples
    G.roots()          dict: component-root triple -> component size
    G.component(root)  set of solution triples in one component

Optional SciPy usage:

    G.eig(root)       second adjacency eigenpair for one component
"""

from __future__ import annotations

import ctypes as _ct
import os
import platform
from pathlib import Path
from typing import Dict, Iterator, List, Set, Tuple

Triple = Tuple[int, int, int]

__all__ = ["MarkoffGraph", "Triple"]
__version__ = "0.0.6"

_MAX_PRIME_BY_NODE_COUNT = 46340

_ERROR_MESSAGES = {
    -1: "invalid argument",
    -2: "too many solutions for the internal node bound",
    -3: "allocation failure inside libmarkoff",
    -4: "internal error: Vieta neighbor was not found among solutions",
    -5: "prime is too large for uint16 coordinates or uint32 node count",
    -6: "modulus is not prime; this builder uses field arithmetic",
    -7: "internal error: too many solutions share a fixed coordinate pair",
    -8: "internal error while building component data or CSR arrays",
}


class _CNode(_ct.Structure):
    _fields_ = [
        ("x", _ct.c_uint16),
        ("y", _ct.c_uint16),
        ("z", _ct.c_uint16),
        ("pad", _ct.c_uint16),
        ("root", _ct.c_uint32),
        ("neighbor0", _ct.c_uint32),
        ("neighbor1", _ct.c_uint32),
        ("neighbor2", _ct.c_uint32),
        ("eigenvector", _ct.c_double),
    ]


class _CComponent(_ct.Structure):
    pass


_CComponent._fields_ = [
    ("root", _ct.POINTER(_CNode)),
    ("root_index", _ct.c_uint32),
    ("size", _ct.c_uint32),
    ("eigenvalue", _ct.c_double),
]


class _CCSR(_ct.Structure):
    _fields_ = [
        ("size", _ct.c_uint32),
        ("nnz", _ct.c_uint32),
        ("root_index", _ct.c_uint32),
        ("data", _ct.POINTER(_ct.c_double)),
        ("indices", _ct.POINTER(_ct.c_int)),
        ("indptr", _ct.POINTER(_ct.c_int)),
        ("nodes", _ct.POINTER(_ct.c_uint32)),
    ]


def _library_names() -> List[str]:
    system = platform.system().lower()
    if system == "windows":
        return ["libmarkoff.dll", "markoff.dll"]
    if system == "darwin":
        return ["libmarkoff.dylib"]
    return ["libmarkoff.so"]


def _load_library(path: str | os.PathLike[str] | None = None) -> _ct.CDLL:
    if path is not None:
        return _ct.CDLL(str(path))

    env_path = os.environ.get("MARKOFF_LIB")
    if env_path:
        return _ct.CDLL(env_path)

    package_dir = Path(__file__).resolve().parent
    errors = []
    for name in _library_names():
        candidate = package_dir / name
        if candidate.exists():
            try:
                return _ct.CDLL(str(candidate))
            except OSError as exc:
                errors.append(f"{candidate}: {exc}")

    searched = ", ".join(str(package_dir / name) for name in _library_names())
    detail = "\n".join(errors)
    raise OSError(
        "Could not load the markoff_graph native library. "
        "Install a platform wheel, build the package locally, set MARKOFF_LIB, "
        f"or pass lib_path. Searched: {searched}"
        + (f"\nLoad errors:\n{detail}" if detail else "")
    )


def _configure_library(lib: _ct.CDLL) -> _ct.CDLL:
    lib.markoff_build.argtypes = [
        _ct.c_int,
        _ct.c_int,
        _ct.c_int,
        _ct.c_int,
        _ct.c_int,
        _ct.POINTER(_ct.c_void_p),
    ]
    lib.markoff_build.restype = _ct.c_int

    lib.markoff_free.argtypes = [_ct.c_void_p]
    lib.markoff_free.restype = None

    lib.markoff_node_count.argtypes = [_ct.c_void_p]
    lib.markoff_node_count.restype = _ct.c_uint32

    lib.markoff_component_count.argtypes = [_ct.c_void_p]
    lib.markoff_component_count.restype = _ct.c_uint32

    lib.markoff_nodes.argtypes = [_ct.c_void_p]
    lib.markoff_nodes.restype = _ct.POINTER(_CNode)

    lib.markoff_components.argtypes = [_ct.c_void_p]
    lib.markoff_components.restype = _ct.POINTER(_CComponent)

    lib.markoff_component_csr.argtypes = [
        _ct.c_void_p,
        _ct.c_uint32,
        _ct.POINTER(_ct.c_void_p),
    ]
    lib.markoff_component_csr.restype = _ct.c_int

    lib.markoff_csr_free.argtypes = [_ct.c_void_p]
    lib.markoff_csr_free.restype = None

    return lib


def _is_prime(n: int) -> bool:
    n = int(n)
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def _normalize(vertex: Triple, p: int) -> Triple:
    x, y, z = vertex
    return int(x) % p, int(y) % p, int(z) % p


def _triple(node: _CNode) -> Triple:
    return int(node.x), int(node.y), int(node.z)


class MarkoffGraph:
    """Connected components of a Markoff-type surface over a prime field.

    The surface is

        x^2 + y^2 + z^2 = x*y*z + A*x + B*y + C*z + D  mod p.

    The modulus must be prime. Coordinates are stored as uint16_t and node
    counts as uint32_t, so this builder requires 2*p*p <= 2^32 - 1.  Thus p
    must be at most 46340; since p is prime, the largest possible p is 46337.
    """

    __slots__ = (
        "_A",
        "_B",
        "_C",
        "_D",
        "_prime",
        "_lib",
        "_handle",
        "_nodes",
        "_components",
        "_node_count",
        "_component_count",
        "_root_sizes",
        "_root_indices",
    )

    def __init__(
        self,
        A: int,
        B: int,
        C: int,
        D: int,
        prime: int,
        *,
        lib_path: str | os.PathLike[str] | None = None,
    ) -> None:
        p = int(prime)
        if p <= 1:
            raise ValueError("prime must be > 1")
        if p > _MAX_PRIME_BY_NODE_COUNT:
            raise ValueError(
                f"prime must be <= {_MAX_PRIME_BY_NODE_COUNT} because 2*p*p must fit in uint32; got {p}"
            )
        if not _is_prime(p):
            raise ValueError(f"modulus must be prime for this field-based builder; got {p}")

        self._A = int(A)
        self._B = int(B)
        self._C = int(C)
        self._D = int(D)
        self._prime = p
        self._lib = _configure_library(_load_library(lib_path))
        self._handle = _ct.c_void_p()
        self._nodes = None
        self._components = None
        self._node_count = 0
        self._component_count = 0
        self._root_sizes: Dict[Triple, int] = {}
        self._root_indices: Dict[Triple, int] = {}

        status = self._lib.markoff_build(
            self._A,
            self._B,
            self._C,
            self._D,
            self._prime,
            _ct.byref(self._handle),
        )
        if status < 0:
            message = _ERROR_MESSAGES.get(status, "unknown error")
            raise RuntimeError(f"markoff_build failed with code {status}: {message}")
        if not self._handle.value:
            raise RuntimeError("markoff_build succeeded but returned a null graph handle")

        self._node_count = int(self._lib.markoff_node_count(self._handle))
        self._component_count = int(self._lib.markoff_component_count(self._handle))
        self._nodes = self._lib.markoff_nodes(self._handle)
        self._components = self._lib.markoff_components(self._handle)
        if self._node_count and not self._nodes:
            self.close()
            raise RuntimeError("markoff_nodes returned null for a nonempty graph")
        if self._component_count and not self._components:
            self.close()
            raise RuntimeError("markoff_components returned null for a nonempty graph")

        self._index_roots()

    def _require_open(self) -> None:
        if self._handle is None or not self._handle.value:
            raise RuntimeError("MarkoffGraph is closed")

    def _node_at(self, i: int) -> _CNode:
        self._require_open()
        assert self._nodes is not None
        return self._nodes[int(i)]

    def _component_at(self, i: int) -> _CComponent:
        self._require_open()
        assert self._components is not None
        return self._components[int(i)]

    def _index_roots(self) -> None:
        root_sizes: Dict[Triple, int] = {}
        root_indices: Dict[Triple, int] = {}

        for i in range(self._component_count):
            component = self._component_at(i)
            root_index = int(component.root_index)
            root = _triple(self._node_at(root_index))
            root_sizes[root] = int(component.size)
            root_indices[root] = root_index

        self._root_sizes = root_sizes
        self._root_indices = root_indices

    def close(self) -> None:
        """Release the C-owned graph memory."""
        if self._handle is not None and self._handle.value:
            self._lib.markoff_free(self._handle)
            self._handle = _ct.c_void_p()
            self._nodes = None
            self._components = None

    def nodes(self) -> Iterator[Triple]:
        """Iterate over all solution triples."""
        self._require_open()
        for i in range(self._node_count):
            yield _triple(self._node_at(i))

    def roots(self) -> Dict[Triple, int]:
        """Return component-root triples with their component sizes."""
        self._require_open()
        return dict(self._root_sizes)

    def component(self, root: Triple) -> Set[Triple]:
        """Return the component whose root is the given root triple."""
        self._require_open()
        root = _normalize(root, self._prime)
        if root not in self._root_indices:
            raise KeyError(f"not a component root modulo {self._prime}: {root!r}")

        root_index = self._root_indices[root]
        return {
            _triple(self._node_at(i))
            for i in range(self._node_count)
            if int(self._node_at(i).root) == root_index
        }


    def _component_root_index(self, root: Triple) -> int:
        self._require_open()
        root = _normalize(root, self._prime)
        if root not in self._root_indices:
            raise KeyError(f"not a component root modulo {self._prime}: {root!r}")
        return self._root_indices[root]

    def _component_csr_arrays(self, root: Triple):
        """Return NumPy arrays for the component CSR adjacency matrix.

        This private helper asks C to build the SciPy-style CSR buffers:
        double* data, int* indices, int* indptr, and a uint32_t* local-to-global
        node map. The buffers are copied into NumPy arrays before the C CSR
        object is freed, so the returned arrays are safe to keep.
        """
        self._require_open()
        try:
            import numpy as np
        except ImportError as exc:
            raise ImportError(
                'CSR export requires numpy. Install it with: python -m pip install "markoff-graph[eig]"'
            ) from exc

        root_index = self._component_root_index(root)
        csr_handle = _ct.c_void_p()
        status = self._lib.markoff_component_csr(
            self._handle,
            _ct.c_uint32(root_index),
            _ct.byref(csr_handle),
        )
        if status < 0:
            message = _ERROR_MESSAGES.get(status, "unknown error")
            raise RuntimeError(f"markoff_component_csr failed with code {status}: {message}")
        if not csr_handle.value:
            raise RuntimeError("markoff_component_csr succeeded but returned a null CSR handle")

        try:
            csr = _ct.cast(csr_handle, _ct.POINTER(_CCSR)).contents
            size = int(csr.size)
            nnz = int(csr.nnz)
            data = np.ctypeslib.as_array(csr.data, shape=(nnz,)).copy()
            indices = np.ctypeslib.as_array(csr.indices, shape=(nnz,)).copy()
            indptr = np.ctypeslib.as_array(csr.indptr, shape=(size + 1,)).copy()
            global_indices = np.ctypeslib.as_array(csr.nodes, shape=(size,)).copy()
            triples = tuple(_triple(self._node_at(int(i))) for i in global_indices)
            return data, indices, indptr, global_indices, triples
        finally:
            self._lib.markoff_csr_free(csr_handle)

    def _store_eig_values(self, root: Triple, eigenvalue: float, global_indices, vector) -> None:
        """Store a SciPy-computed eigenpair back into the C-owned graph object."""
        self._require_open()
        root_index = self._component_root_index(root)
        assert self._components is not None
        assert self._nodes is not None

        for ci in range(self._component_count):
            component = self._components[ci]
            if int(component.root_index) == root_index:
                component.eigenvalue = float(eigenvalue)
                break

        for gi, value in zip(global_indices, vector):
            self._nodes[int(gi)].eigenvector = float(value)

    def eig(self, root: Triple):
        """Return ``(eigenvalue, eigenvector_dict)`` for one component.

        This uses the optional SciPy backend. Install it with:

            python -m pip install "markoff-graph[eig]"
        """
        from .eig import eig

        return eig(self, root)

    def __enter__(self) -> "MarkoffGraph":
        self._require_open()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self) -> str:
        status = "closed" if self._handle is None or not self._handle.value else "open"
        return (
            f"MarkoffGraph(nodes={self._node_count}, components={self._component_count}, "
            f"prime={self._prime}, {status})"
        )
