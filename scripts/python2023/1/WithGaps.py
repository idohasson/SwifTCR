import csv
import itertools
import sys
import itertools
from collections import defaultdict

import numpy as np
import pandas as pd
from tempfile import mkstemp

def write_cluster(cdr3_buffer, result_filename):
    results = []
    for cid, cdr3aa in cdr3_buffer.items():
        results += [[aa, cid] for aa in cdr3aa]
    pd.DataFrame(results, columns=['CDR3', 'cluster_id']).to_csv(result_filename, mode="a", index=False)
    # with open(output_path, 'a', newline='') as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerows([[cdr3_seq, cluster_id] for cdr3_seq in cdr3_list])
    #     csvfile.flush()

if __name__ == "__main__":

    WRITE_BUFFER = 100

    file_path = sys.argv[1] if len(sys.argv) > 1 else "listCDR3.csv"
    result_path = sys.argv[2] if len(sys.argv) > 2 else "output.csv"
    row_limit = int(sys.argv[3]) if len(sys.argv) > 3 else None

    # Read file
    print("Reading file")
    aa_cdr3s = pd.read_csv(file_path, nrows=row_limit, header=0, usecols=["aaSeqCDR3"], squeeze=True)
    aa_cdr3s = aa_cdr3s.unique()

    # create hashes
    hash_seqs = []
    cdr3_i = []
    for i, aa in enumerate(aa_cdr3s):
        # find all minus-one-length combinations of the CDR3
        hash_frag = set("".join(a) for a in itertools.combinations(aa, len(aa) - 1))
        # hash each combination
        hash_seqs += [hash(aa)] + list(hash(frag) for frag in hash_frag)
        cdr3_i += [i] * (len(hash_frag) + 1)


    hash_seqs, cdr3_i = np.array(hash_seqs), np.array(cdr3_i)

    sorted_i = hash_seqs.argsort()

    hash_seqs, cdr3_i = hash_seqs[sorted_i], cdr3_i[sorted_i]

    # find cliques
    clique_l = []
    c_indices = []
    last = hash_seqs[0]
    for h, i in zip(hash_seqs, cdr3_i):
        if last == h:  # same hash for same string
            c_indices.append(i)
        elif len(c_indices) > 1:
            clique_l.append(list(c_indices))
            c_indices = [i]
        else:
            c_indices = [i]
        last = h

    # # find all connecting sequences
    unique_i = list(set(id for ids in clique_l for id in ids))
    unique_i = {j: i for i, j in enumerate(unique_i)}

    # create edges
    edge_list = defaultdict(list)
    for clique in clique_l:
        aa_i = [unique_i[ci] for ci in clique]
        edges = list(itertools.permutations(aa_i, 2))
        for u, v in edges:
            edge_list[u].append(v)

    # BFS
    visited = [False] * len(unique_i)
    queue = []
    w_buffer = defaultdict(list)
    cid = 0

    while not all(visited):
        visit_next = next((i for i, x in enumerate(visited) if not x))
        queue.append(visit_next)

        while queue:
            v = queue.pop(0)
            if not visited[v]:
                visited[v] = True
                queue.extend([u for u in edge_list[v] if not visited[u]])
                w_buffer[cid].append(aa_cdr3s[v])

            if len(w_buffer) > WRITE_BUFFER:
                write_cluster(w_buffer, result_path)
                w_buffer.clear()
        cid += 1

    write_cluster(w_buffer, result_path)
    w_buffer.clear()


    pass




