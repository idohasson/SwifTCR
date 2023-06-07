
# -------------------- network_func --------------------

import random
import string
import time


def seq(length):
    aa = string.ascii_uppercase[:20]
    return ''.join(random.choice(aa) for _ in range(length))


def seq_list(num, min_len=10, max_len=17):
    return [seq(random.randint(min_len, max_len)) for _ in range(num)]


# def randint_list_unique(a, b, n):
#     return set(random.randint(a, b, n))


# generate a graph of cliques from a list of sets so edges are colored based on the number of set they are in
def generate_graph_cliques(sets):
    G = nx.Graph()
    for clique in sets:
        for i, j in combinations(clique, 2):
            # check if the edge is already in the graph
            if not G.has_edge(i, j):
                G.add_edge(i, j, shared=0)
            else:
                G[i][j]['shared'] += 1
    return G


# sum clique edges to get the number of sets they are in
def sum_clique_edges(clique_network):
    shared = nx.get_edge_attributes(clique_network, 'shared').values()
    return sum(shared)


# draw the graph so edges are colored based on the set they are in
def draw_graph(draw_network):
    colors = nx.get_edge_attributes(draw_network, 'shared').values()
    nx.draw(draw_network, with_labels=True, node_color='#A0CBE2', edge_color=colors, cmap=cm.Blues, vmin=0, vmax=2)
    plt.show()


# draw dense graph with colored edges based on the number of sets they are in
def draw_dense_graph(dense_network):
    shared = nx.get_edge_attributes(dense_network, 'shared').values()
    nx.draw(dense_network, with_labels=True, node_color='#A0CBE2', edge_color=shared, cmap=cm.Blues, vmin=0, vmax=2)
    plt.show()


def timit(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f'{func.__name__} took {end - start:.3f} seconds')
        return result

    return wrapper


@timit
def test(f, s_list=None, times=10):
    for _ in range(times):
        if isinstance(s_list, int):
            x = [f(seq, 1) for seq in seq_list(s_list)]
        else:
            random.shuffle(s_list)
            x = [f(seq, 1) for seq in s_list]


# --------------------  --------------------

from itertools import combinations
from numpy import random
import networkx as nx
from matplotlib import pyplot as plt
from matplotlib import cm


def randint_list_unique(a, b, n):
    return set(random.randint(a, b, n))


# generate a graph of cliques from a list of sets so edges are colored based on the number of set they are in
def generate_graph_cliques(sets):
    G = nx.Graph()
    for clique in sets:
        for i, j in combinations(clique, 2):
            # check if the edge is already in the graph
            if not G.has_edge(i, j):
                G.add_edge(i, j, shared=0)
            else:
                G[i][j]['shared'] += 1
    return G


# sum clique edges to get the number of sets they are in
def sum_clique_edges(clique_network):
    shared = nx.get_edge_attributes(clique_network, 'shared').values()
    return sum(shared)


# draw the graph so edges are colored based on the set they are in
def draw_graph(draw_network):
    colors = nx.get_edge_attributes(draw_network, 'shared').values()
    nx.draw(draw_network, with_labels=True, node_color='#A0CBE2', edge_color=colors, cmap=cm.Blues, vmin=0, vmax=2)
    plt.show()


# draw dense graph with colored edges based on the number of sets they are in
def draw_dense_graph(dense_network):
    shared = nx.get_edge_attributes(dense_network, 'shared').values()
    nx.draw(dense_network, with_labels=True, node_color='#A0CBE2', edge_color=shared, cmap=cm.Blues, vmin=0, vmax=2)
    plt.show()
