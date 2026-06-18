# markoff-graph
build A/B/C/D markoff graph modulo p
```python
from markoff_graph import markoff

G = markoff(4, 4, -2, -4, 31)

G.nodes()
G.edges()
G.roots()
G.component()
G.neighbors((x, y, z))
```
