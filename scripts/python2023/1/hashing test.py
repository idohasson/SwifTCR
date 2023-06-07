import gc
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

import itertools
from IPython.display import clear_output


def get_seq_combinations(full_sequences):
    return set(itertools.combinations(full_sequences, len(full_sequences) - 1))


def get_hashes(combinations_sets):
    merged_sets = itertools.chain(*combinations_sets)
    hash_sets = map(hash, merged_sets)
    return hash_sets


def read_file(file_path, row_number=1000):
    n = 0
    pp = int(row_number / 100)

    with open(file_path, "r") as csvfile:
        datareader = csv.reader(csvfile)
        next(datareader)  # yield the header row
        for row in datareader:
            yield row

            n += 1
            if not n % pp:
                clear_output(wait=True)
                print(np.round(100 * n / row_number, 2), "%", end="", flush=True)
            if n == row_number:
                clear_output(wait=True)
                return


def get_hash_comb(sequence_reader):
    sil = []
    hl = []
    n = 0
    for row, seq in sequence_reader:
        comb = get_seq_combinations(seq)
        h_combs = list(map(hash, comb))
        hl += h_combs

        combs_n = len(h_combs)
        sil += [seq] * combs_n
        n += combs_n

    return hl, sil


def sorted_indices(to_sort):
    return np.argsort(np.array(to_sort))


def get_cliques(hl, sorted_i):
    clique_s = []
    clique = set()
    for si1, si2 in itertools.zip_longest(sorted_i[:-1], sorted_i[1:]):
        if hl[si1] == hl[si2]:
            clique.add(si1)
            clique.add(si2)
        elif clique:
            clique_s.append(clique.copy())
            clique.clear()
    return clique_s


def get_edges(from_clique, seq_i):
    em = defaultdict(list)
    for cs in from_clique:
        si_cs = [seq_i[oi] for oi in cs]
        e = list(itertools.permutations(si_cs, 2))
        for u, v in e:
            em[u].append(v)
    return em


def find_clusters(edges):
    visited = dict.fromkeys(edges.keys(), False)
    queue = [list(edges.keys())[0]]
    cid = 0
    w_cbuffer = []
    while queue:
        v = queue.pop(0)
        if not visited[v]:
            visited[v] = True
            queue += [u for u in edges[v] if not visited[u]]
            w_cbuffer.append((v, cid))

        if not queue:
            if not all(visited.values()):
                visit_next = next(i for i, x in visited.items() if not x)
                queue.append(visit_next)
            cid += 1
    return w_cbuffer


def write_clusters(w_path, w_clusters):
    with open(w_path, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(["seq_id", "cluster_id"])
        # write the data
        for cid, w_cluster in enumerate(w_clusters):
            for wc in w_cluster:
                writer.writerow([wc, cid])


if __name__ == '__main__':
    import csv
    import time

    FILE_PATH = sys.argv[1]
    LIMIT_READ = int(sys.argv[2])

    # seq_i_list = []
    # hash_list = []
    # n = 0
    t0 = time.time()

    seq_reader = read_file(FILE_PATH, LIMIT_READ)

    hash_list, seq_i_list = get_hash_comb(seq_reader)

    # for row, seq in read_file(FILE_PATH, LIMIT_READ):
    #     # if not isinstance(seq, str):
    #     #     continue
    #     # if len(seq) < 5:
    #     #     continue
    #     # process row
    #     comb = get_seq_combinations(seq)
    #     h_combs = list(map(hash, comb))
    #     hash_list += h_combs
    #
    #     combs_n = len(h_combs)
    #     seq_i_list += [int(row)] * combs_n
    #     n += combs_n

    # sorted_indices = np.argsort(np.array(hash_list))
    sorted_indices = sorted_indices(hash_list)
    # t1 = time.time()
    # print(np.round(t1 - t0, 3), "seconds for", LIMIT_READ, "sequences.")

    # find cliques
    # clique_sets = []
    clique_sets = get_cliques(hash_list, sorted_indices)

    # clique = set()
    # for si1, si2 in itertools.zip_longest(sorted_indices[:-1], sorted_indices[1:]):
    #     if hash_list[si1] == hash_list[si2]:
    #         clique.add(si1)
    #         clique.add(si2)
    #     elif clique:
    #         clique_sets.append(clique.copy())
    #         clique.clear()

    # BFS
    ## create adges

    edge_map = get_edges(clique_sets, seq_i_list)
    # edge_map = defaultdict(list)
    # for cs in clique_sets:
    #     si_cs = [seq_i_list[oi] for oi in cs]
    #     edges = list(itertools.permutations(si_cs, 2))
    #     for u, v in edges:
    #         edge_map[u].append(v)

    clusters = find_clusters(edge_map)

    # visited = dict.fromkeys(edge_map.keys(), False)
    # queue = [list(edge_map.keys())[0]]
    # cid = 0
    # w_buffer = []
    # while queue:
    #     v = queue.pop(0)
    #     if not visited[v]:
    #         visited[v] = True
    #         queue += [u for u in edge_map[v] if not visited[u]]
    #         w_buffer.append((v, cid))
    #
    #     if not queue:
    #         if not all(visited.values()):
    #             visit_next = next(i for i, x in visited.items() if not x)
    #             queue.append(visit_next)
    #         cid += 1

    write_clusters('test_write2.csv', clusters)
    #
    # with open('test_write.csv', 'w', encoding='UTF8') as f:
    #     writer = csv.writer(f)
    #     # write the header
    #     writer.writerow(["seq_id", "cluster_id"])
    #     # write the data
    #     writer.writerows(clusters)
