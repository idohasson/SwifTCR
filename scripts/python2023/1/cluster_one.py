import numpy as np
from itertools import groupby
from operator import itemgetter


def hamming_cluster(sequences):
    """
    Take a list of sequences and cluster them by hamming distance 1.
    @param sequences - the list of sequences to cluster
    @return the clusters
    """

    hashes = np.array([(hash(s[:i] + s[i + 1:]) ^ i, si)
                       for si, s in enumerate(sequences)
                       for i in range(len(s))])

    hashes = hashes[hashes[:, 0].argsort(), :]
    hashes = np.split(hashes[:, 1], np.unique(hashes[:, 0], return_index=True)[1][1:])

    return [[sequences[i] for i in c] for c in filter(lambda g: g.shape[0] > 1, hashes)]


def levenshtein_cluster(sequences):
    """
    Given a list of sequences, cluster them by the Levenshtein distance 1.
    @param sequences - the list of sequences to cluster           
    @return The list of clusters           
    """
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
