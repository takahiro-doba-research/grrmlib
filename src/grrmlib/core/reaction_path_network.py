import copy
from typing import Callable, Hashable, Self

import networkx as nx

from .molecule import Molecule
from .molecules import Molecules
from .grouped_network import GroupedNetwork


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
    
    def add_eqs(self, mols: Molecules) -> Self:
        nodes = [(mol.name, {"EQ": mol}) for mol in mols.values()]
        self.add_nodes_from(nodes)
        return self
    
    def add_pts(self, mols: Molecules) -> Self:
        edges = [(*mol.connection, mol.name, {"PT": mol}) for mol in mols.values()]
        self.add_edges_from(edges)
        return self
    
    def copy(self) -> Self:
        return copy.deepcopy(self)
    
    def remove_eq(self, name: Hashable) -> Self:
        if name in self:
            self.remove_node(name)
        return self
    
    def get_pts(self, u: Hashable, v: Hashable) -> Molecules:
        edge_data = self.get_edge_data(u, v)
        
        if edge_data is None:
            return Molecules()
        
        return Molecules(
            {k: d["PT"].copy() for k, d in edge_data.items()}
        )
    
    def cluster(
        self,
        predicate: Callable[[Molecule, Molecule], bool]
    ) -> GroupedNetwork:
        pass