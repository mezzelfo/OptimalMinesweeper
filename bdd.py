import matplotlib.pyplot as plt
import networkx as nx
from itertools import product

def getMyBDD(n_vars, target):
    def yolo(thisnode):
        level, partial = thisnode
        node_sonTrue = (level+1, partial+1)
        node_sonFalse = (level+1, partial)
        # Gestisco l'ultimo livello. Ci dovrebbero essere solo due nodi al fondo
        if level + 1 == n_vars:
            if partial == target:
                G.add_edge(thisnode, 1, value = 0)
                G.add_edge(thisnode, 0, value = 1)
            elif partial+1 == target:
                G.add_edge(thisnode, 0, value = 0)
                G.add_edge(thisnode, 1, value = 1)
            else:
                raise RuntimeError("Nope")
        else:
            # Gestisto tutti gli edge 1
            if partial + 1 > target:
                    G.add_edge(thisnode, 0, value = 1)
            else:
                G.add_edge(thisnode, node_sonTrue, value = 1)
                yolo(node_sonTrue)

            # Gestisto tutti gli edge 0
            if (n_vars-level-1)+partial < target:
                G.add_edge(thisnode, 0, value = 0)
            else:
                G.add_edge(thisnode, node_sonFalse, value = 0)
                yolo(node_sonFalse)

    G = nx.DiGraph()
    yolo((0,0))
    return G

for N in range(2,10):
    for T in range(N):
        G = getMyBDD(N,T)

        for path in nx.all_simple_edge_paths(G, source = (0,0), target = 0):
            vals = [G.get_edge_data(u,v)['value'] for (u,v) in path]
            assert sum(vals) != T

        for path in nx.all_simple_edge_paths(G, source = (0,0), target = 1):
            vals = [G.get_edge_data(u,v)['value'] for (u,v) in path]
            try:
                assert sum(vals) == T
            except Exception:
                print(path,vals,sum(vals))
        print(N,T,G.number_of_nodes(),G.number_of_edges())
