import time
from collections import defaultdict
from itertools import chain, combinations, count, filterfalse, repeat

from more_itertools import prepend, split_when

from SwifTCR.hash_token import timeit


# from memory_profiler import profile

def timit(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f'{func.__name__} took {end - start:.3f} seconds')
        return result

    return wrapper



def levenshtein(sequence_list):
    '''
    Given a list of sequences, return a list of sets of sequences that are all identical except for one character.

    Args:
    sequence_list: a list of sequences to be compared
    Returns:
    A list of sets of sequences that are all identical except for one character.
    '''

    def hash_subsequences(seq):
        sub = prepend(tuple(seq), combinations(seq, len(seq) - 1))
        return zip(map(hash, sub), count(), repeat(seq))

    def hash_clusters(groups):
        clusters = [(c[0][1], set(map(itemgetter(-1), c)))
                    for c in split_when(groups, lambda s1, s2: s1[1] != s2[1])]
        if clusters[0][0] == 0:
            clusters = [(p, c.union(clusters[0][1])) for p, c in clusters[1:]]
        return filterfalse(lambda g: len(g[1]) == 1, clusters)

    sub_keys = chain(*map(hash_subsequences, sequence_list))
    sorted_keys = sorted(sub_keys, key=lambda k: k[0] + k[1])

    hash_groups = split_when(sorted_keys, lambda s1, s2: s1[0] != s2[0])
    hash_groups = filterfalse(lambda g: len(g) == 1, hash_groups)
    hash_groups = chain(*map(hash_clusters, hash_groups))

    return hash_groups



@timeit
def create_combinations_dict(sequences):
    d = defaultdict(set)
    for x, seq in enumerate(sequences):
        for y, z in zip(combinations(seq, len(seq) - 1), zip(count(), repeat(x))):
            d[y].add(z)
    return d


# ---------------------------- #
def cdr3_group(list_size, sequence_i):
    return np.searchsorted(np.cumsum(list_size), sequence_i)


# ---------------------------- #

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
