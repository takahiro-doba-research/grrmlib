import networkx as nx

from .molecules import PTs


class ReactionPathNetwork:
    
    def __init__(self):
        self._G = nx.MultiGraph()
    
    def __getattr__(self, name):
        return getattr(self._G, name)
    
    def to_networkx(self):
        return self._G
    
    def add_eqs(self, eqs):
        nodes = [(eq.name, {"EQ": eq}) for eq in eqs.values()]
        self._G.add_nodes_from(nodes)
    
    def add_pts(self, pts):
        edges = [(*pt.connection, pt.name, {"PT": pt}) for pt in pts.values()]
        self._G.add_edges_from(edges)
    
    def remove_eq(self, name):
        if name in self._G:
            self._G.remove_node(name)
    
    def get_pts(self, u, v):
        return PTs({n: d["PT"] for n, d in self._G.get_edge_data(u, v).items()})