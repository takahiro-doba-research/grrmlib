import networkx as nx

from .molecules import Molecules


class GroupedNetwork:
    
    def __init__(self):
        self._G = nx.Graph()
    
    def __getattr__(self, name):
        return getattr(self._G, name)
    
    def to_networkx(self):
        return self._G
    
    def add_rpn(self, rpn):
        nodes = []
        for n, eq in rpn.nodes(data="EQ"):
            n_group = eq.group
            for node in nodes:
                if node[0] == n_group:
                    node[1]["EQs"][n] = eq
                    break
            else:
                nodes.append((n_group, {"EQs": Molecules({n: eq})}))
        self._G.add_nodes_from(nodes)
        
        edges = []
        for u, v, k, pt in rpn.edges(keys=True, data="PT"):
            u_group = rpn.nodes[u]["EQ"].group
            v_group = rpn.nodes[v]["EQ"].group
            for edge in edges:
                if set(edge[:2]) == {u_group, v_group}:
                    edge[2]["PTs"][k] = pt
                    break
            else:
                edges.append((u_group, v_group, {"PTs": Molecules({k: pt})}))
        self._G.add_edges_from(edges)
    
    def get_eqs(self, n):
        return self._G.nodes[n]["EQs"]
        
    def get_pts(self, u, v):
        return self._G.edges[u, v]["PTs"]