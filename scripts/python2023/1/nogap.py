import itertools
import sys
import time
import numpy as np
import pandas as pd
import csv
from collections import defaultdict


def create_hash(seq):
    # find all minus-one-length combinations of the CDR3
    seq_frag = set("".join(a) for a in itertools.combinations(seq, len(seq) - 1))
    # hash each combination
    frag_hash = [hash(frag) for frag in seq_frag]
    return frag_hash


def read_cdr3_file(fp, limit=5000):
    df = pd.read_csv(fp, nrows=limit, header=0, index_col=0)
    cdr3s = [row[0] for index, row in df.iterrows()]
    return cdr3s


def write_cluster(cdr3_list, cluster_id, output_path):
    with open(output_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows([[cdr3_seq, cluster_id] for cdr3_seq in cdr3_list])


def find_cliques(s_hashes):
    clique_l = []
    cq = []
    last = s_hashes[0]
    for i, h in enumerate(s_hashes):
        if last == h:  # same hash for same string
            cq.append(i)
        elif len(cq) > 1:
            clique_l.append(cq)
            cq = [i]
        else:
            cq = [i]
        last = h
    return clique_l


def create_edges(vertexes):
    edge_list = []
    for vertex in vertexes:
        edge_list.extend(list(itertools.permutations(vertex, 2)))
    return edge_list


class Graph:

    def __init__(self, vertices):
        self.V = vertices  # No. of vertices
        self.graph = defaultdict(list)  # default dictionary to store graph

    # function to add an edge to graph
    def addEdge(self, u, v):
        self.graph[u].append(v)

    # A function used by DFS
    def DFSUtil(self, v, visited):
        # Mark the current node as visited and print it
        visited[v] = True

        # Recur for all the vertices adjacent to this vertex
        for i in self.graph[v]:
            if visited[i] == False:
                self.DFSUtil(i, visited)
        return visited

    def fillOrder(self, v, visited, stack):
        # Mark the current node as visited
        visited[v] = True
        # Recur for all the vertices adjacent to this vertex
        for i in self.graph[v]:
            if visited[i] == False:
                self.fillOrder(i, visited, stack)
        stack = stack.append(v)

    # Function that returns reverse (or transpose) of this graph
    def getTranspose(self):
        g = Graph(self.V)

        # Recur for all the vertices adjacent to this vertex
        for i in self.graph:
            for j in self.graph[i]:
                g.addEdge(j, i)
        return g

    # The main function that finds and prints all strongly
    # connected components
    def printSCCs(self):

        stack = []
        # Mark all the vertices as not visited (For first DFS)
        visited = [False] * (self.V)
        # Fill vertices in stack according to their finishing
        # times
        for i in range(self.V):
            if visited[i] == False:
                self.fillOrder(i, visited, stack)

        # Create a reversed graph
        gr = self.getTranspose()

        # Mark all the vertices as not visited (For second DFS)
        visited = [False] * (self.V)

        # Now process all vertices in order defined by Stack
        components_collection = []
        while stack:
            i = stack.pop()
            if visited[i] == False:
                new_visit = [False] * (self.V)
                component = gr.DFSUtil(i, new_visit)
                ci = [i for i, c in enumerate(component) if c]
                for i in ci:
                    visited[i] = True
                components_collection.append(ci)

        return components_collection


def get_cc(seq_list):
    hash_seqs = []
    seq_ids = []
    for i, seq in enumerate(seq_list):
        hash_frags = create_hash(seq)
        hash_seqs.extend(hash_frags)
        seq_ids.extend([i] * len(hash_frags))

    hash_seqs = np.array(hash_seqs)
    sorted_i = hash_seqs.argsort()
    sorted_seq = hash_seqs[sorted_i]

    clique_list = find_cliques(sorted_seq)
    seq_ids = np.array(seq_ids)[sorted_i]

    clique_ids = [seq_ids[ids] for ids in clique_list]

    clique_list = [[seq_list[sid] for sid in s_ids] for s_ids in clique_ids]

    unique_s = []
    for clique in clique_list:
        unique_s.extend([cl for cl in clique])
    unique_s = list(set(unique_s))

    clique_ids = [[unique_s.index(sc) for sc in clique] for clique in clique_list]

    el = create_edges(clique_ids)

    g = Graph(len(unique_s))
    for v, u in el:
        g.addEdge(v, u)
    scc = g.printSCCs()
    # if scc:
    #     all_components[len_key] = [[unique_seq[index] for index in cc] for cc in scc]
    scc = [[unique_s[cci] for cci in cc] for cc in scc]
    return scc


if __name__ == "__main__":
    start_time = time.time()  # TIMESTAMP
    time_file_name = "time.csv"
    WRITE_BUFFER = 100
    file_path = sys.argv[1] if len(sys.argv) > 1 else "listCDR3.csv"
    result_file_name = sys.argv[2] if len(sys.argv) > 2 else "output.csv"
    row_limit = int(sys.argv[3]) if len(sys.argv) > 3 else None

    # Read file
    print("Reading file")

    seq_list = read_cdr3_file(file_path, row_limit)

    read_file_time = time.time()  # TIMESTAMP

    cdr3_dict = defaultdict(list)
    for seq in seq_list:
        if len(seq) > 5:
            cdr3_dict[len(seq)].append(seq)

    connected_component = [get_cc(cdr3s) for cdr3s in cdr3_dict.values()]

    # print("Hashing")
    # hashes = [(len(seq), create_hash(seq)) for seq in seq_list if len(seq) > 5]
    # hash_dict = defaultdict(list)
    # cdr3_dict = defaultdict(list)
    # for seq_l, hash_t in hashes:
    #     h, s = [ht[0] for ht in hash_t], [ht[1] for ht in hash_t]
    #     hash_dict[seq_l].extend(h)
    #     cdr3_dict[seq_l].extend(s)
    #
    # sequence_hash_time = time.time()
    #
    # print("Sorting")
    #
    # # sort each hash list separately:
    # for key, hash_list in hash_dict.items():
    #     hash_arr = np.array(hash_list)
    #     si = np.array(hash_list).argsort()
    #     hash_dict[key] = list(hash_arr[si])
    #     cdr3_dict[key] = [cdr3_dict[key][i] for i in si]
    #
    #
    # # Python implementation of Kosaraju's algorithm to print all SCCs
    #
    # # This class represents a directed graph using adjacency list representation
    #
    # def find_cliques(sorted_hashes, sequences):
    #     clique_list = []
    #     clique = []
    #     last = sorted_hashes[0]
    #     for i, h in enumerate(sorted_hashes):
    #         if last == h:  # same hash for same string
    #             clique.append(i)
    #         elif len(clique) > 1:
    #             clique_list.append([sequences[ci] for ci in clique])
    #             clique = [i]
    #         else:
    #             clique = [i]
    #         last = h
    #     return clique_list
    #
    #
    # def create_edges(vertexes):
    #     edge_list = []
    #     for vertex in vertexes:
    #         edge_list.extend(list(itertools.permutations(vertex, 2)))
    #     return edge_list
    #
    #
    # all_components = {}
    # for len_key in hash_dict.keys():
    #     # qi = find_cliques(hash_dict[len_key], cdr3_dict[len_key])
    #     clique_list = find_cliques(hash_dict[len_key], cdr3_dict[len_key])
    #
    #     # clique_list, s = find_cliques(hash_dict[len_key])
    #     # clique_list = find_cliques(hash_dict[len_key])  # TODO - sort by size
    #     # convert aa sequences to int
    #     unique_seq = []
    #     for c in clique_list:
    #         unique_seq.extend(c)
    #     unique_seq = list(set(unique_seq))
    #
    #     seq_to_indices = [[unique_seq.index(seq) for seq in cc] for cc in clique_list]
    #     el = create_edges(seq_to_indices)
    #
    #     g = Graph(len(unique_seq))
    #     for v, u in el:
    #         g.addEdge(v, u)
    #     scc = g.printSCCs()
    #     if scc:
    #         all_components[len_key] = [[unique_seq[index] for index in cc] for cc in scc]

    all_clusters = []
    for ac in connected_component:
        all_clusters.extend(ac)

    for cluster_id, cluster in enumerate(all_clusters):
        write_cluster(cluster, cluster_id, result_file_name)
