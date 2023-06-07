import os
import sys
import time

import numpy as np
import pandas as pd
import itertools
from collections import defaultdict


def write_times(file_path, timestamps, limit_read):
    delta = np.array(list(times.values())[1:]) - np.array(list(times.values())[:-1])
    delta = np.round(delta, 3)
    results = dict(zip(list(times.keys())[1:], delta))
    results["total"] = np.round(time.time() - timestamps["start"], 3)

    if not os.path.exists(file_path):
        pd.DataFrame(results, index=[limit_read]).to_csv(file_path, mode='w', sep='\t')
    else:
        pd.DataFrame(results, index=[limit_read]).to_csv(file_path, mode='a', sep='\t', header=False)


def write_cluster(cdr3_buffer, directory, filename):
    if not os.path.exists(directory):
        os.makedirs(directory)
    o_path = os.path.join(directory, filename)
    cdr3_buffer = [[seq, cluster_id] for cluster_id, cdr3_seq in cdr3_buffer.items() for seq in cdr3_seq]
    pd.DataFrame(cdr3_buffer, columns=['CDR3', 'cluster_id']).to_csv(o_path, mode="a", index=False)

def read_file(file_path):
    print("Reading file")
    return pd.read_csv(file_path, nrows=LIMIT_READ, header=0, usecols=["aaSeqCDR3"], squeeze=True)

def get_hashes(cdr3_list):
    hash_s = []
    cdr3_index = []
    for i, aa in enumerate(cdr3_list):
        # find all minus-one-length combinations of the CDR3
        hash_frag = set("".join(a) for a in itertools.combinations(aa, len(aa) - 1))
        # hash each combination
        hash_s += [hash(aa)] + list(hash(frag) for frag in hash_frag)
        cdr3_index += [i] * (len(hash_frag) + 1)
    return hash_s, cdr3_index


def find_cliques(hash_s, index):
    c_indices = []
    last = hash_s[0]
    edge_l = defaultdict(list)
    for h, i in zip(hash_s, index):
        if last == h:  # same hash for same string
            c_indices.append(i)
        elif len(c_indices) > 1:
            edges = list(itertools.permutations(aa_cdr3s[c_indices], 2))
            for u, v in edges:
                edge_l[u].append(v)
            c_indices = [i]
        else:
            c_indices = [i]
        last = h
    return edge_l

if __name__ == "__main__":
    # assert len(sys.argv) > 2, "First argument: Input file is needed"
    # assert os.path.isfile(sys.argv[1]) and os.access(sys.argv[1], os.R_OK), "First argument: Could not find a CSV file."
    # input_path = sys.argv[1]
    # output_path = sys.argv[2]
    # result_filename = sys.argv[3]
    # LIMIT_READ = int(sys.argv[4]) if len(sys.argv) > 4 else None
    input_path = "listCDR3.csv"
    output_path = "results"
    result_filename = "test2.csv"
    LIMIT_READ = 1000


    WRITE_BUFFER = 2**13
    # timing
    times = {"start": time.time()}


    aa_cdr3s = read_file(input_path)
    # pd.read_csv(input_path, nrows=LIMIT_READ, header=0, usecols=["aaSeqCDR3"], squeeze=True)



    times["reading"] = time.time()
    # create hashes
    print("Hashing")
    hash_seqs, cdr3_index = get_hashes(aa_cdr3s)

    times["hash"] = time.time()
    print("Sorting")

    hash_seqs, cdr3_index = np.array(hash_seqs), np.array(cdr3_index)
    sorted_i = hash_seqs.argsort()
    hash_seqs, cdr3_index = hash_seqs[sorted_i], cdr3_index[sorted_i]

    times["sort"] = time.time()

    # find cliques

    edge_list = find_cliques(hash_seqs, cdr3_index)

    # create edges
    times["edges"] = time.time()

    # BFS
    print("Clustering")
    visited = {aa: False for aa in edge_list.keys()}
    queue = [list(visited.keys())[0]]
    cid = 0
    w_buffer = defaultdict(list)
    while queue:
        v = queue.pop(0)
        if not visited[v]:
            visited[v] = True
            queue.extend([u for u in edge_list[v] if not visited[u]])

            if sys.getsizeof(w_buffer) + sys.getsizeof(v) > WRITE_BUFFER:
                write_cluster(w_buffer, output_path, result_filename)
                w_buffer.clear()

            w_buffer[cid].append(v)

        if not queue:
            if not all(visited.values()):
                visit_next = next(i for i, x in visited.items() if not x)
                queue.append(visit_next)
            cid += 1

    write_cluster(w_buffer, output_path, result_filename)
    w_buffer.clear()

    times["search"] = time.time()

    write_times("running_time.csv", times, LIMIT_READ)


