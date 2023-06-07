import pandas as pd
import csv
from itertools import combinations, chain, groupby
from operator import itemgetter
from collections import Counter
from more_itertools import flatten, sliding_window


def cdr3_read(file_path, min_len=8, max_len=20, field="aaSeqCDR3", sep=","):
    df = pd.read_csv(file_path, index_col=None, header=0, usecols=[field], sep=sep, dtype=str,
                     na_filter=False).drop_duplicates()
    return {CDR3_length: group['aaSeqCDR3'].tolist() for CDR3_length, group in df.groupby(df.aaSeqCDR3.str.len()) if
            min_len <= CDR3_length <= max_len}


def edge_to_string(cluster_edges, from_i):
    rows = []
    for i, edges in enumerate(cluster_edges, from_i):
        for edge in edges:
            if isinstance(edge[-1], tuple):
                if edge[2] == "g":
                    rows.append([edge[0], edge[1], str(i), edge[2], ",".join(str(pos) for pos in edge[3])])
                elif edge[2] == "t":
                    rows.append([edge[0], edge[1], str(i), edge[2], "-".join(str(pos) for pos in edge[3])])
            else:
                rows.append([edge[0], edge[1], str(i), edge[2], str(edge[3])])
    return rows


def write_clusters(records, cluster_id, write_to="output.tsv"):
    write_type = "a"
    if cluster_id == 0:
        write_type = "w"
    rows = edge_to_string(records, cluster_id)
    with open(write_to, write_type, newline='') as outcsv:
        writer = csv.writer(outcsv, delimiter='\t')
        if cluster_id == 0:
            writer.writerow(["CDR3_1", "CDR3_2", "cluster_id", "edit_type", "edit_position"])
        writer.writerows(rows)


def hash_tokens(seq, seq_len, seq_n):
    hash_list = []
    for i, s in enumerate(seq):
        if i < seq_n:
            for comb in combinations(s, seq_len - 1):
                hash_list.append(hash(comb))
        else:
            hash_list.append(hash(tuple(s)))
    return hash_list

##############  path to file with ##############
# assumes file has header
# assumes CDR3 in column name 'aaSeqCDR3' (for different column name use argument: field="column_name")
path = r"C:\Users\hassoni4\PycharmProjects\testing\data\standard_MiXCR.txt"
# path = "/work/datasets/file_sampling/size_100000/sample_100000_0.csv"

##############  for comma separated values (csv) ##############
# cdr3_sequences = cdr3_read(path, min_len=5, max_len=25)

##############  for tab separated values  ##############
cdr3_sequences = cdr3_read(path, 5, 25, sep="\t")

tokens_n = map(len, cdr3_sequences.values())
tokens = map(tuple, flatten(cdr3_sequences.values()))
search_order = list(sliding_window(chain([0], cdr3_sequences.keys()), 2))

cid = 0
for l_shorter, l in search_order:
    print(l)
    cdr3 = []
    cdr3 += cdr3_sequences[l]
    n = len(cdr3)

    if l_shorter:
        cdr3 += cdr3_sequences[l_shorter]

    hashes = hash_tokens(cdr3, l, n)
    si = sorted(range(len(hashes)), key=lambda i: hashes[i])

    edge_list, t_singleton = [], set()
    for k, g in groupby(enumerate(si), lambda ix: hashes[ix[1]]):
        cluster, indels = [], {}
        for i, j in combinations(map(itemgetter(1), g), 2):

            if int(i / l) == int(j / l):
                continue

            if max(i, j) / l >= n:
                ii, jj = sorted([i, j])
                indels.setdefault((cdr3[int(ii / l)], cdr3[jj - (n - 1) * l]), []).append(l - ii % l - 1)

            elif i % l == j % l:
                cluster.append((cdr3[int(i / l)], cdr3[int(j / l)], "s", l - i % l - 1))

            elif abs(i % l - j % l) == 1:
                ii, jj = i - i % l + j % l, j - j % l + i % l
                if hashes[ii] == hashes[jj] and hashes[ii] != hashes[i] and hashes[j] != hashes[jj]:
                    ii, jj = sorted([i, j], key=lambda ij: ij % l, reverse=True)
                    cluster.append((cdr3[int(i / l)], cdr3[int(j / l)], "t", (l - ii % l - 1, l - jj % l - 1)))

        if indels:
            for aa_pair, edit_pos in indels.items():
                if len(edit_pos) > 1:
                    cluster.append((aa_pair[0], aa_pair[1], "i", tuple(sorted(edit_pos))))
                else:
                    cluster.append((aa_pair[0], aa_pair[1], "i", edit_pos[0]))

        if len(cluster) == 1 and cluster[0][2] == "t":
            if cluster[0] not in t_singleton:
                t_singleton.add(cluster.pop(0))

        if cluster:
            edge_list.append(cluster)

    for ts, tsn in Counter(t_singleton).items():
        if tsn == 2:
            edge_list.append([ts])

    write_clusters(edge_list, cid)
    cid += len(edge_list)

    edge_list.clear()
    t_singleton.clear()
