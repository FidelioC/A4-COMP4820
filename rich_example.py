import graphviz

dot = graphviz.Digraph(format="svg")
prev = graph[0]
for node in graph[1:]:
    dot.node(prev)
    dot.node(node)
    dot.edge(prev, node)
    prev = node
