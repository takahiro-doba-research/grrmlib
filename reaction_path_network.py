import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from pyvis.network import Network

from .geometry import get_adj_matrix


class ReactionPathNetwork:
    
    def __init__(self, eq_list, pt_list):
        self.eq_list = eq_list
        self.pt_list = pt_list
        
        G = nx.MultiGraph()
        nodes = [(eq.name, eq.__dict__) for eq in eq_list]
        G.add_nodes_from(nodes)
        edges = [(*pt.connection, pt.__dict__) for pt in pt_list]
        G.add_edges_from(edges)
        if "EQ??" in G:
            G.remove_node("EQ??")
        self.rpn = G
    
    def to_html(self, path):
        # Creat another graph H for visualization.
        H = self.rpn.copy()
        
        # Set node attributes.
        nx.set_node_attributes(H, {n: {"id": n} for n in H.nodes()})
        nx.set_node_attributes(H, {n: {"label": n} for n in H.nodes()})
        nx.set_node_attributes(H, {n: {"title": n} for n in H.nodes()})
        nx.set_node_attributes(H, 5, "size")
        nx.set_node_attributes(H, "dot", "shape")
        
        energy_eq0 = H.nodes["EQ0"]["scfenergy"]
        nx.set_node_attributes(
            H, {n: {"scfenergy_adjusted": (data["scfenergy"] - energy_eq0) * 2625.5} for n, data in H.nodes(data=True)}
        )
        
        def energy2color(energy, cmap, norm):
            color = cmap(norm(energy))
            color = tuple(c * s for c, s in zip(color, (255, 255, 255, 1)))
            color = f"rgb{color}"
            return color
        
        cmap = plt.cm.cividis
        norm = Normalize(-100, 0)
        nx.set_node_attributes(
            H, {n: {"color": energy2color(data["scfenergy_adjusted"], cmap, norm)} for n, data in H.nodes(data=True)}
        )
        
        groups = {}  # dict of dict
        
        for node, data in H.nodes(data=True):
            arr_adj = get_adj_matrix(data["symbols"], data["atomcoords"])
            for key, value in groups.items():
                if np.all(value["adj_matrix"] == arr_adj):
                    value["group_nth"] += 1
                    H.nodes[node]["group"] = key
                    H.nodes[node]["group_nth"] = value["group_nth"]
                    break
            else:
                key = f"G{len(groups)}"
                H.nodes[node]["group"] = key
                H.nodes[node]["group_nth"] = 0
                groups[key] = ({"group_nth": 0, "adj_matrix":arr_adj})

        n_groups = len(groups)
        nx.set_node_attributes(
            H, {
                n: {
                    "x": 10 * n_groups * (2 + data["group_nth"]) * np.cos(2 * np.pi * int(data["group"][1:]) / n_groups),
                    "y": - 10 * n_groups * (2 + data["group_nth"]) * np.sin(2 * np.pi * int(data["group"][1:]) / n_groups),
                }
                for n, data in H.nodes(data=True)
            }
        )
        
        # Set edge attributes.
        nx.set_edge_attributes(
            H, {(u, v, k): {"label": data["name"]} for u, v, k, data in H.edges(keys=True, data=True)}
        )
        nx.set_edge_attributes(
            H, {(u, v, k): {"title": data["name"]} for u, v, k, data in H.edges(keys=True, data=True)}
        )
        nx.set_edge_attributes(
            H, {(u, v, k): {"smooth": {"type": "curvedCW", "roundness": 0.1 + 0.1 * k}} for u, v, k in H.edges(keys=True)}
        )
        nx.set_edge_attributes(
            H, {(u, v, k): {"arrows": {"to": {"enabled": False}}} for u, v, k in H.edges(keys=True)}
        )
        nx.set_edge_attributes(
            H, {
                (u, v, k): {"scfenergy_adjusted": (data["scfenergy"] - energy_eq0) * 2625.5}
                for u, v, k, data in H.edges(keys=True, data=True)
            }
        )
        nx.set_edge_attributes(
            H, {
                (u, v, k): {"color": energy2color(data["scfenergy_adjusted"], cmap, norm)}
                for u, v, k, data in H.edges(keys=True, data=True)
            }
        )

        # Creat pyvis network.
        net = Network(directed=True, notebook=True, height="100%", width="100%")
        
        for node, data in H.nodes(data=True):
            data = {k: v for k, v in data.items() if k in ["id", "label", "title", "size", "shape", "x", "y", "color"]}
            net.add_node(node, **data)
            
        for source, target, key, data in H.edges(keys=True, data=True):
            data = {k: v for k, v in data.items() if k in ["label", "title", "smooth", "arrows", "color"]}
            net.add_edge(source, target, **data)
        
        for i in range(n_groups):
            net.add_node(
                f"G{i}",
                **{
                    "id": f"G{i}",
                    "label": f"G{i}",
                    "title": f"G{i}",
                    "shape": "circle",
                    "color": {"background": "rgba(0,0,0,0)", "border": "black"},
                    "x": 10 * n_groups * np.cos(2 * np.pi * i / n_groups),
                    "y": - 10 * n_groups * np.sin(2 * np.pi * i / n_groups),
                }
            )
        
        net.toggle_physics(False)
        net.set_edge_smooth("dynamic")
        net.show(str(path))
    
    def to_html2(self, path):
        # Creat another graph H for visualization.
        H = self.rpn.copy()
        
        # Set node attributes.
        nx.set_node_attributes(H, {n: {"id": n} for n in H.nodes()})
        nx.set_node_attributes(H, {n: {"label": n} for n in H.nodes()})
        nx.set_node_attributes(H, {n: {"title": n} for n in H.nodes()})
        nx.set_node_attributes(H, 5, "size")
        nx.set_node_attributes(H, "dot", "shape")
        
        energy_eq0 = H.nodes["EQ0"]["scfenergy"]
        nx.set_node_attributes(
            H,{n: {"scfenergy_adjusted": (data["scfenergy"] - energy_eq0) * 2625.5} for n, data in H.nodes(data=True)}
        )
        
        def energy2color(energy, cmap, norm):
            color = cmap(norm(energy))
            color = tuple(c * s for c, s in zip(color, (255, 255, 255, 1)))
            color = f"rgb{color}"
            return color
        
        cmap = plt.cm.cividis
        norm = Normalize(-100, 0)
        nx.set_node_attributes(
            H, {node: {"color": energy2color(data["scfenergy_adjusted"], cmap, norm)} for node, data in H.nodes(data=True)}
        )
        
        n_nodes = len(H)
        nx.set_node_attributes(
            H, {
                n: {
                    "x": 10 * n_nodes * np.cos(2 * np.pi * i / n_nodes),
                    "y": - 10 * n_nodes * np.sin(2 * np.pi * i / n_nodes),
                }
                for i, n in enumerate(H.nodes())
            }
        )
        
        # Set edge attributes.
        nx.set_edge_attributes(
            H, {(u, v, k): {"label": data["name"]} for u, v, k, data in H.edges(keys=True, data=True)}
        )
        nx.set_edge_attributes(
            H, {(u, v, k): {"title": data["name"]} for u, v, k, data in H.edges(keys=True, data=True)}
        )
        nx.set_edge_attributes(
            H, {(u, v, k): {"smooth": {"type": "curvedCW", "roundness": 0.1 + 0.1 * k}} for u, v, k in H.edges(keys=True)}
        )
        nx.set_edge_attributes(
            H, {(u, v, k): {"arrows": {"to": {"enabled": False}}} for u, v, k in H.edges(keys=True)}
        )
        nx.set_edge_attributes(
            H, {
                (u, v, k): {"scfenergy_adjusted": (data["scfenergy"] - energy_eq0) * 2625.5}
                for u, v, k, data in H.edges(keys=True, data=True)
            }
        )
        nx.set_edge_attributes(
            H, {
                (u, v, k): {"color": energy2color(data["scfenergy_adjusted"], cmap, norm)}
                for u, v, k, data in H.edges(keys=True, data=True)
            }
        )
        
        # Creat pyvis network.
        net = Network(directed=True, notebook=True, height="100%", width="100%")
        
        for node, data in H.nodes(data=True):
            data = {k: v for k, v in data.items() if k in ["id", "label", "title", "size", "shape", "x", "y", "color"]}
            net.add_node(node, **data)
            
        for source, target, key, data in H.edges(keys=True, data=True):
            data = {k: v for k, v in data.items() if k in ["label", "title", "smooth", "arrows", "color"]}
            net.add_edge(source, target, **data)
        
        net.toggle_physics(False)
        net.set_edge_smooth("dynamic")
        net.show(str(path))