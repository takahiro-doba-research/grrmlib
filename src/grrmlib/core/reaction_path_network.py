from typing import Callable, Hashable, Self

import networkx as nx

from .molecule import Molecule
from .molecules import Molecules
from .grouped_reaction_path_network import GroupedReactionPathNetwork


class ReactionPathNetwork(nx.MultiGraph):
    
    def __init__(
        self,
        *,
        eqs: Molecules | None = None,
        pts: Molecules | None = None,
    ) -> None:
        super().__init__()
        
        if eqs is not None:
            self.add_eqs(eqs)
        
        if pts is not None:
            self.add_pts(pts)
    
    def add_eqs(self, eqs: Molecules) -> Self:
        nodes = [(eq.name, {"EQ": eq}) for eq in eqs.values()]
        self.add_nodes_from(nodes)
        return self
    
    def add_pts(self, pts: Molecules) -> Self:
        edges = [(*pt.connection, pt.name, {"PT": pt}) for pt in pts.values()]
        self.add_edges_from(edges)
        return self
    
    def remove_eq(self, name: Hashable) -> Self:
        if name in self:
            self.remove_node(name)
        return self
    
    def pts_between(self, u: Hashable, v: Hashable) -> Molecules:
        edge_data = self.get_edge_data(u, v)
        
        if edge_data is None:
            return Molecules()
        
        return Molecules(
            {k: d["PT"].copy() for k, d in edge_data.items()}
        )
    
    def cluster(
        self,
        predicate: Callable[[Molecule, Molecule], bool]
    ) -> GroupedReactionPathNetwork:
        grpn = GroupedReactionPathNetwork()
        
        # Add grouped EQs
        eqs = Molecules(dict(self.nodes(data="EQ")))
        geqs = eqs.cluster(predicate)
        nodes = [(g, {"EQs": eqs}) for g, eqs in geqs.items()]
        grpn.add_nodes_from(nodes)
        
        # Add PTs
        eq_to_group = {
            name: group
            for group, eqs in geqs.items()
            for name in eqs
        }
        edges: dict[frozenset | tuple, Molecules] = {}
        
        for u, v, name, pt in self.edges(keys=True, data="PT"):
            u_group = eq_to_group[u]
            v_group = eq_to_group[v]
            
            if u_group == v_group:
                key = (u_group, v_group)
            else:
                key = frozenset((u_group, v_group))
            
            if key not in edges:
                edges[key] = Molecules()
            
            edges[key][name] = pt.copy()
        
        edges = [(*key, {"PTs": pts}) for key, pts in edges.items()]
        grpn.add_edges_from(edges)
        
        return grpn