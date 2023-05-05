import csv
import sys
import itertools
import time

import numpy as np
from collections import defaultdict
from memory_profiler import profile
from IPython.display import clear_output


def read_file(file_path, col_i=1):
    # N = 0
    # pp = int(row_number / 100)
    with open(file_path, "r") as csvfile:
        datareader = csv.reader(csvfile)
        next(datareader)  # yield the header row
        for row in datareader:
            aa_seq = row[col_i]
            if isinstance(aa_seq, str):
                aa_seq_len = len(aa_seq)
                if aa_seq_len > 5:
                    yield aa_seq_len, aa_seq

            # N += 1
            # if not N % pp:
            #     clear_output(wait=True)
            #     print("Reading...", np.round(100 * N / row_number, 2), "%", end="", flush=True)
            # if N == row_number:
            #     return


def get_comb(sequence_reader):

    def hash_combinations(full_seq, seq_l):
        sub_seq_comb = list(itertools.combinations(full_seq, seq_l - 1))
        subs_set = ["".join(sub_seq) for sub_seq in sub_seq_comb]
        encoded = [sub.encode("utf-8") for sub in subs_set]
        hash_subs_set = [hash(e_seq) for e_seq in encoded]
        hash_subs_set = list(set(hash_subs_set))
        return hash_subs_set

    full_seq = defaultdict(list)
    for seq_l, seq in sequence_reader:
        full_seq[seq_l].append(seq)

    for seq_len, f_seqs in full_seq.items():
        comb_sets = [hash_combinations(seq, seq_len) for seq in f_seqs]
        hashses = np.array([c for c in itertools.chain.from_iterable(comb_sets)])
        hashses_i = np.argsort(hashses)
        sets_s = np.array([full_seq[seq_len][i] for i, hs in enumerate(comb_sets) for _ in range(len(hs))])
        yield sets_s[hashses_i], hashses[hashses_i]


def sorted_indices(to_sort):
    return np.argsort(np.array(to_sort))


def get_cliques(seq_list, ordered_i):
    clique = set()
    hash_pairs = itertools.zip_longest(ordered_i[:-1], ordered_i[1:])
    seq_pairs = itertools.zip_longest(seq_list[:-1], seq_list[1:])

    for (si1, si2), (h1, h2) in zip(seq_pairs, hash_pairs):
        if h1 == h2:
            clique.add(si1)
            clique.add(si2)
        elif clique:
            yield clique.copy()
            clique.clear()


def get_edges(from_clique):
    em = defaultdict(set)
    for clique in from_clique:
        for v in clique:
            em[v].update(clique)
            em[v].discard(v)
    return em


# def connected_component(edges):
#     visited = dict.fromkeys(edges.keys(), False)
#     queue = [list(edges.keys())[0]]
#     w_cbuffer = [[]]
#     while queue:
#         v = queue.pop()
#         if not visited[v]:
#             visited[v] = True
#             queue += [u for u in edges[v] if not visited[u]]
#             w_cbuffer[-1].append(v)
#
#         if not queue:
#             if not all(visited.values()):
#                 for i, x in visited.items():
#                     if not x:
#                         queue.append(i)
#                         break
#                 for j in w_cbuffer[-1]:
#                     del visited[j]
#                     del edges[j]  # check if needed
#                 w_cbuffer.append([])
#     return w_cbuffer

def connected_component(edges):
    left = set(edges.keys())
    queue = set()
    queue.add(left.pop())
    cluster = set()
    all_clusters = []
    while queue:
        v = queue.pop()
        if not cluster:
            cluster.add(v)

        neighbours = edges[v].difference(cluster)
        if neighbours:
            cluster.update(neighbours)
            queue.update(neighbours)

        if not queue:
            # write
            all_clusters.append(cluster.copy())
            left.difference_update(cluster)

            if left:
                next_v = left.pop()
                queue.add(next_v)
                cluster.clear()
    return all_clusters


def find_clusters(sorted_hashes):
    for st, oi in sorted_hashes:
        clique_s = get_cliques(st, oi)
        edge_map = get_edges(clique_s)
        if edge_map:
            yield connected_component(edge_map)


def write_clusters(w_path, clusters_gen):
    n = 0
    for w_clusters in clusters_gen:
        if not n:
            mode = 'w'
            write_header = True
        else:
            mode = 'a'
            write_header = False
        wc = []

        for cluster_id, cluster_set in enumerate(w_clusters, n):
            wc += [(aa, cluster_id) for aa in cluster_set]
            n += 1

        with open(w_path, mode, encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            # write the header
            if write_header:
                writer.writerow(["CDR3", "cluster_id"])
            # write the data
            writer.writerows(wc)
            for w in wc:
                print(w)


# @profile
def main(input_path, output_path):
    seq_reader = read_file(input_path)

    hash_combinations = get_comb(seq_reader)

    cluster_list = find_clusters(hash_combinations)

    write_clusters(output_path, cluster_list)


if __name__ == '__main__':
    i_path = sys.argv[1]
    o_path = sys.argv[2]
    main(i_path, o_path)
