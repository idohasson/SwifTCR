import numpy as np
import pandas as pd
from itertools import combinations, count, groupby, repeat, starmap
from operator import getitem, itemgetter
from more_itertools import unique_justseen
from memprof import memprof



def hamming_cluster(sequences):
    hashes = np.array([(hash(s[:i] + s[i + 1:]) ^ i, si)
                       for si, s in enumerate(sequences)
                       for i in range(len(s))])

    hashes = hashes[hashes[:, 0].argsort(), :]
    hashes = np.split(hashes[:, 1], np.unique(hashes[:, 0], return_index=True)[1][1:])

    return [[sequences[i] for i in c] for c in filter(lambda g: g.shape[0] > 1, hashes)]


def levenshtein_cluster(sequences):
    hashes = np.array([(hash(s[:i] + s[i + 1:]) ^ (len(s) - i), si)
                       for si, s in enumerate(sequences) for i in range(len(s) + 1)])

    hashes = hashes[hashes[:, 0].argsort(), :]

    unique_h, index = np.unique(hashes[:, 0], return_index=True)
    index = [index[i] for i in range(1, len(unique_h)) if
             abs(unique_h[i] - unique_h[i - 1]) > 10]  # CHANGE TO MAX SEQUENCE LENGTH INSTEAD OF 10

    clusters = []
    for grouped_i in np.split(hashes, index, axis=0):

        cluster = [set(sequences[i] for i in map(itemgetter(1), g))
                   for _, g in groupby(grouped_i.tolist(), key=itemgetter(0))]

        if len(set(map(len, set.union(*cluster)))) == 2:
            gap = min(set.union(*cluster), key=len)
            for i in range(len(cluster)): cluster[i].add(gap)

        clusters += [c for c in cluster if len(c) > 1]

    return clusters


def hash_subsequence(sequences, skip):
    hash_token = lambda s: hash(s[:skip] + s[skip + 1:])
    return np.fromiter(map(hash_token, sequences), dtype=np.int64)


def hash_cluster(sequences, skip_index, with_i):
    hashes = hash_subsequence(sequences[with_i:], skip_index)

    clusters = groupby(hashes.argsort(), key=lambda i: getitem(hashes, i))

    clusters = starmap(lambda h, g: [with_i + i for i in unique_justseen(g, key=lambda i: sequences[with_i + i])],
                       clusters)

    clusters = filter(lambda c: len(c) > 1, clusters)

    return clusters


def distance_one_hamming(input_cdr3, cdr3_index=False):
    input_cdr3.sort(key=len)
    pos, full_s = 0, 0

    for to_i, max_l in unique_justseen(enumerate(map(len, input_cdr3)), key=itemgetter(1)):

        while pos < max_l:

            clusters = hash_cluster(input_cdr3, pos, to_i)

            if not cdr3_index:
                clusters = map(lambda c: [input_cdr3[i] for i in c], clusters)

            yield from zip(clusters, repeat(pos))

            full_s, pos = to_i, pos + 1


def edge_list(clustered_sequences, write_to_file=None):
    edge_list_df = pd.DataFrame()
    cid = count()
    for cluster, pos in clustered_sequences:

        clustered_edge_list_df = pd.DataFrame(combinations(cluster, 2), columns=['CDR3_1', 'CDR3_2'])
        clustered_edge_list_df['cluster_id'] = next(cid)
        clustered_edge_list_df['edit_position'] = pos
        edge_list_df = pd.concat([edge_list_df, clustered_edge_list_df], ignore_index=True)

    if write_to_file is not None:
        edge_list_df.to_csv(write_to_file, index=False)

    return edge_list_df
