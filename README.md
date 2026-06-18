# markoff-graph

Small Python package for finite-field Markoff-type graphs over prime fields:

```text
x^2 + y^2 + z^2 = x*y*z + A*x + B*y + C*z + D  mod p
```

## Install

```bash
python -m pip install markoff-graph
```

Optional features:

```bash
python -m pip install "markoff-graph[nx]"    # NetworkX export
python -m pip install "markoff-graph[eig]"   # SciPy eigenpairs
python -m pip install "markoff-graph[all]"   # all optional features
```

## Usage

```python
from markoff_graph import MarkoffGraph

G = MarkoffGraph(4, 4, -2, -4, 31)

print(G.roots())                 # root triple -> component size
print(len(G.component((0, 0, 19))))
print(list(G.nodes())[:5])
```

## NetworkX export

```python
from markoff_graph import MarkoffGraph
from markoff_graph.nx_export import to_multigraph

G = MarkoffGraph(4, 4, -2, -4, 31)
H = to_multigraph(G)
```

## Eigenpair for a component

```python
from markoff_graph import MarkoffGraph

G = MarkoffGraph(2, 1, 0, -5, 41)
value, vector = G.eig((0, 0, 6))

print(value)          # second adjacency eigenvalue
print(vector[(0, 0, 6)])
```

