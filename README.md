# markoff-graph
build A/B/C/D markoff graph modulo p
```python
>>> from markoff_graph import markoff
>>> G = markoff(4, 4, -2, -4, 31)
>>> print(G)
MarkoffGraph(nodes=962, components=2, prime=31)
>>> print(G.roots())
>>> for (x, y, z) in G.roots():
...     print((x, y, z), G.neighbors((x, y, z)))
...     print('component size', len(G.component((x, y, z))))
...     
(0, 0, 19) ((4, 0, 19), (0, 4, 19), (0, 0, 10))
component size 450
(0, 1, 30) ((3, 1, 30), (0, 3, 30), (0, 1, 30))
component size 512
```
