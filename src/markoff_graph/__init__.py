"""Graphs of Markoff-type surfaces over prime fields.

Usage:

    from markoff_graph import markoff

    G = markoff(4, 4, -2, -4, 31)
    print(G.roots())

The graph object exposes:

    G.nodes()                 iterator of solution triples (x, y, z)
    G.edges()                 iterator of directed edges ((x,y,z), (x',y',z'))
    G.roots()                 dict: component-root triple -> component size
    G.component()             dict: component-root triple -> set of node triples
    G.neighbors((x, y, z))    3-tuple of Vieta-neighbor triples
"""

from __future__ import annotations

import ctypes as _ct
import os
import platform
from array import array
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

Triple = Tuple[int, int, int]

__all__ = ["Triple", "markoff"]
__version__ = "0.1.0"

_MAX_PRIME_U32 = 1621
_U32_ARRAY_CODE = "I"

_ERROR_MESSAGES = {
    -1: "invalid input pointer, prime <= 1, or capacity too large",
    -2: "too many solutions for the allocated uint32 buffers",
    -3: "allocation failure inside libmarkoff",
    -4: "internal error: Vieta neighbor was not found among solutions",
    -5: f"prime > {_MAX_PRIME_U32}; uint32 vertex encoding may overflow",
    -6: "modulus is not prime; this builder uses field arithmetic",
    -7: "internal error: too many solutions share a fixed coordinate pair",
}


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
    lib.build_graph.argtypes = [
        _ct.c_int, _ct.c_int, _ct.c_int, _ct.c_int,
        _ct.c_int,
        _ct.c_uint32,
        _ct.POINTER(_ct.c_int),
        _ct.POINTER(_ct.c_uint32),
        _ct.POINTER(_ct.c_uint32),
    ]
    lib.build_graph.restype = _ct.c_int
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


def _u32_buffer(capacity: int) -> array:
    if array(_U32_ARRAY_CODE).itemsize != _ct.sizeof(_ct.c_uint32):
        raise RuntimeError("array('I') is not 32-bit on this platform")
    return array(_U32_ARRAY_CODE, [0]) * int(capacity)


def _as_u32_ptr(buf: array):
    return (_ct.c_uint32 * len(buf)).from_buffer(buf)


def _default_capacity(p: int) -> int:
    return 8 if p == 2 else 2 * int(p) * int(p)


def _decode_vertex(V: int, p: int) -> Triple:
    z, rem = divmod(int(V), p * p)
    y, x = divmod(rem, p)
    return int(x), int(y), int(z)


def _normalize(vertex: Triple, p: int) -> Triple:
    x, y, z = vertex
    return int(x) % p, int(y) % p, int(z) % p


class _MarkoffGraph:
    __slots__ = (
        "_A", "_B", "_C", "_D", "_prime",
        "_nodes_u32", "_roots_u32", "_root_sizes", "_component_count",
    )

    def __init__(
        self,
        *,
        A: int,
        B: int,
        C: int,
        D: int,
        prime: int,
        nodes_u32: array,
        roots_u32: array,
        component_count: int,
    ) -> None:
        self._A = int(A)
        self._B = int(B)
        self._C = int(C)
        self._D = int(D)
        self._prime = int(prime)
        self._nodes_u32 = nodes_u32
        self._roots_u32 = roots_u32
        self._component_count = int(component_count)

        root_sizes: Dict[Triple, int] = {}
        for r in self._roots_u32:
            root = _decode_vertex(int(r), self._prime)
            root_sizes[root] = root_sizes.get(root, 0) + 1
        self._root_sizes = root_sizes

    def nodes(self) -> Iterator[Triple]:
        p = self._prime
        for V in self._nodes_u32:
            yield _decode_vertex(int(V), p)

    def edges(self) -> Iterator[Tuple[Triple, Triple]]:
        for v in self.nodes():
            for w in self.neighbors(v):
                yield v, w

    def roots(self) -> Dict[Triple, int]:
        return dict(self._root_sizes)

    def component(self) -> Dict[Triple, Set[Triple]]:
        p = self._prime
        components: Dict[Triple, Set[Triple]] = {}
        for V, R in zip(self._nodes_u32, self._roots_u32):
            root = _decode_vertex(int(R), p)
            node = _decode_vertex(int(V), p)
            components.setdefault(root, set()).add(node)
        return components

    def neighbors(self, vertex: Triple) -> Tuple[Triple, Triple, Triple]:
        p = self._prime
        x, y, z = _normalize(vertex, p)
        if not self._is_solution((x, y, z)):
            raise KeyError(f"not a solution vertex modulo {p}: {vertex!r}")

        sx = ((y * z + self._A - x) % p, y, z)
        sy = (x, (x * z + self._B - y) % p, z)
        sz = (x, y, (x * y + self._C - z) % p)
        return sx, sy, sz

    def _is_solution(self, vertex: Triple) -> bool:
        p = self._prime
        x, y, z = vertex
        lhs = (x*x + y*y + z*z) % p
        rhs = (x*y*z + self._A*x + self._B*y + self._C*z + self._D) % p
        return lhs == rhs

    def __repr__(self) -> str:
        return (
            f"MarkoffGraph(nodes={len(self._nodes_u32)}, "
            f"components={self._component_count}, prime={self._prime})"
        )


def markoff(
    A: int,
    B: int,
    C: int,
    D: int,
    prime: int,
    *,
    lib_path: str | os.PathLike[str] | None = None,
    capacity: Optional[int] = None,
) -> _MarkoffGraph:
    """Build the Vieta graph modulo a prime p.

    The graph is formed from solutions to

        x^2 + y^2 + z^2 = x*y*z + A*x + B*y + C*z + D  mod p.

    The modulus must be prime and at most 1621.
    """
    p = int(prime)
    if p <= 1:
        raise ValueError("prime must be > 1")
    if p > _MAX_PRIME_U32:
        raise ValueError(f"prime must be <= {_MAX_PRIME_U32}; got {p}")
    if not _is_prime(p):
        raise ValueError(f"modulus must be prime for this field-based builder; got {p}")

    cap = int(_default_capacity(p) if capacity is None else capacity)
    if cap <= 0:
        raise ValueError("capacity must be positive")

    lib = _configure_library(_load_library(lib_path))
    nodes = _u32_buffer(cap)
    roots = _u32_buffer(cap)
    n_out = _ct.c_int(0)

    component_count = lib.build_graph(
        int(A), int(B), int(C), int(D), p, _ct.c_uint32(cap),
        _ct.byref(n_out), _as_u32_ptr(nodes), _as_u32_ptr(roots),
    )
    if component_count < 0:
        message = _ERROR_MESSAGES.get(component_count, "unknown error")
        extra = f"; capacity={cap}, partial_count={n_out.value}" if component_count == -2 else ""
        raise RuntimeError(f"build_graph failed with code {component_count}: {message}{extra}")

    n = int(n_out.value)
    del nodes[n:]
    del roots[n:]

    return _MarkoffGraph(
        A=int(A), B=int(B), C=int(C), D=int(D), prime=p,
        nodes_u32=nodes, roots_u32=roots,
        component_count=int(component_count),
    )
