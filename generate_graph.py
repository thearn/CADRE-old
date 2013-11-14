from CADRE.CADRE_assembly import CADRE
import networkx as nx
import numpy as np
import pprint

assembly = CADRE(150, 100)

graph = assembly._depgraph

defaults = ["derivative_exec_count", "directory", "itername", "exec_count",
            "force_execute", "driver"]

remove = []
for node in graph.nodes_iter():
    for d in defaults:
        if d in node:
            remove.append(node)
            break

for node in remove:
    graph.remove_node(node)

ag = nx.to_agraph(graph)
ag.layout('dot')
ag.draw('design.pdf')
