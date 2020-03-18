"""
   A Module of utility functions to 
   visualize CI Engine related steps  
"""

import networkx as nx
import matplotlib.pyplot as plt
from base import UCIEngine
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def visualize_1e_inter_cat_interactions(CIEng, ax=None, layout='spring', show_edges=False, options=None):
    assert isinstance(CIEng, UCIEngine)
    assert CIEng.Cats['graph_1e'] is not None
    
    graph_1e = CIEng.Cats['graph_1e']
    G = nx.DiGraph(directed=True)
    G.add_nodes_from(list(range(CIEng.NCat)))

    ex_levels = {'REF':[0]}
    nodes = set([0])
    cur_nodes = nodes.copy()
    labels = ['S','D','T','Q','5','6','7','8','9','10']
    
    label_i = 0
    while len(nodes) != CIEng.NCat:
        next_nodes = set()
        for i in cur_nodes:
            for j in graph_1e[i].keys():
                if j not in nodes:
                    nodes.add(j)
                    next_nodes.add(j)
    
        cur_nodes = next_nodes.copy()
        ex_levels[labels[label_i]] = (list(cur_nodes))
        label_i += 1

    if show_edges: edge_labels = {}
    
    for i in range(CIEng.NCat):
        for j in graph_1e[i].keys():
            if i==j: continue
            G.add_edge(i,j)

            if show_edges and (graph_1e[i][j][0]<graph_1e[i][j][1]):
                edge_labels[(i,j)] = str([graph_1e[i][j][k]+1 for k in (0,1)])

    if ax is None: fig, ax = plt.subplots(figsize=(10,10))

    if options is None:
        options = {
            'node_size': 1000,
            'width': 3,
            'arrowstyle': '-|>',
            'arrowsize': 12,
            'alpha': 1.0,
            'font_size':20
        }

    if layout=='spring':
        pos = nx.spring_layout(G)
    elif layout=='shell':
        pos = nx.shell_layout(G, list(ex_levels.values()))
    else:
        raise Exception()

    color_map = ['C']*CIEng.NCat
    i = 0
    legend_elements = []
    for ex in ex_levels.keys():
        for node_j in ex_levels[ex]:
            color_map[node_j] += str(i)
        
        legend_elements.append(Line2D([0], [0], marker='o', color='w', label=ex,
                               markerfacecolor='C'+str(i), markersize=20))
        i +=1

    # [color_map[node] for node in G.nodes()] if nodes are not in order

    nx.draw(G, pos, ax=ax, arrows=True, node_color=color_map, 
            labels={node:node for node in G.nodes()},**options)
    
    if show_edges:
        nx.draw_networkx_edge_labels(G,pos,ax=ax,edge_labels=edge_labels,
                                     font_color='red',font_size=10)

    ax.legend(handles=legend_elements, loc='best', bbox_to_anchor=(1.1, 1.1),
              labelspacing=1.25,fontsize=15)

    return ax




    
    




