from typing import Hashable, Self

import networkx as nx

from .molecules import Molecules


class GroupedReactionPathNetwork(nx.Graph):
    
    def eqs_of_group(self, group: Hashable) -> Molecules:
        return self.nodes[group]["EQs"].copy()
    
    def group_of_eq(self, eq_name: Hashable) -> Hashable | None:
        for group, eqs in self.nodes(data="EQs"):
            if eq_name in eqs:
                return group
        return None
    
    def pts_between(self, group1: Hashable, group2: Hashable) -> Molecules:
        edge_data = self.get_edge_data(group1, group2)
        
        if edge_data is None:
            return Molecules()
        
        return edge_data.get("PTs", Molecules()).copy()