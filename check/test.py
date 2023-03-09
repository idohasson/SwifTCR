from collections import defaultdict
from itertools import chain, combinations, count, repeat, groupby, starmap
from operator import itemgetter
from typing import List

import numpy as np

# from SwifTCR.hash_token import timeit


# # from memory_profiler import profile
from SwifTCR.all.random_data import rand_rep

# # def timit(func):
# #     def wrapper(*args, **kwargs):
# #         start = time.time()
# #         result = func(*args, **kwargs)
# #         end = time.time()
# #         print(f'{func.__name__} took {end - start:.3f} seconds')
# #         return result
# #
# #     return wrapper
# #
# #
# #
# # def levenshtein(sequence_list):
# #     '''
# #     Given a list of sequences, return a list of sets of sequences that are all identical except for one character.
# #
# #     Args:
# #     sequence_list: a list of sequences to be compared
# #     Returns:
# #     A list of sets of sequences that are all identical except for one character.
# #     '''
# #
# #     def hash_subsequences(seq):
# #         sub = prepend(tuple(seq), combinations(seq, len(seq) - 1))
# #         return zip(map(hash, sub), count(), repeat(seq))
# #
# #     def hash_clusters(groups):
# #         clusters = [(c[0][1], set(map(itemgetter(-1), c)))
# #                     for c in split_when(groups, lambda s1, s2: s1[1] != s2[1])]
# #         if clusters[0][0] == 0:
# #             clusters = [(p, c.union(clusters[0][1])) for p, c in clusters[1:]]
# #         return filterfalse(lambda g: len(g[1]) == 1, clusters)
# #
# #     sub_keys = chain(*map(hash_subsequences, sequence_list))
# #     sorted_keys = sorted(sub_keys, key=lambda k: k[0] + k[1])
# #
# #     hash_groups = split_when(sorted_keys, lambda s1, s2: s1[0] != s2[0])
# #     hash_groups = filterfalse(lambda g: len(g) == 1, hash_groups)
# #     hash_groups = chain(*map(hash_clusters, hash_groups))
# #
# #     return hash_groups
# #
# # def create_combinations_dict(sequences):
# #     d = defaultdict(set)
# #     for x, seq in enumerate(sequences):
# #         for y, z in zip(combinations(seq, len(seq) - 1), zip(count(), repeat(x))):
# #             d[y].add(z)
# #     return d
# #
# # # ---------------------------- #
# # def cdr3_group(list_size, sequence_i):
# #     return np.searchsorted(np.cumsum(list_size), sequence_i)
# # # ---------------------------- #
# # import numpy as np
# # from itertools import groupby
# # from operator import itemgetter
# #


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


def hamming_cluster2(sequences):
    """
    Take a list of sequences and cluster them by hamming distance 1.
    @param sequences - the list of sequences to cluster
    @return the clusters
    """

    hashes = np.array([(hash(s[:i] + s[i + 1:]) + i, si)
                       for si, s in enumerate(sequences)
                       for i in range(len(s))], dtype=[('hash', np.int64), ('pos', np.int64)])
    hashes = np.sort(hashes, order='hash')
    hashes = np.split(hashes, np.unique(hashes['hash'], return_index=True)[1][1:])
    return [[sequences[i] for i in c['pos']] for c in filter(lambda g: g.shape[0] > 1, hashes)]


def hamming_cluster3(sequences):
    """
    Take a list of sequences and cluster them by hamming distance 1.
    @param sequences - the list of sequences to cluster
    @return the clusters
    """

    # def hash_subsequence(seq):
    #     '''Return L-1 length subsequences of the input string.'''
    #     return map(hash, combinations(seq, len(seq) - 1))
    #
    # hashed_tokens = enumerate(map(hash_subsequence, sequences))
    # hashed_tokens = chain(*map(lambda x: zip(x[1], repeat(x[0])), hashed_tokens))
    # # hash tokens function

    def hash_token(idx, seq):
        tokens = combinations(seq, len(seq) - 1)
        tokens = map(hash, tokens)
        tokens = zip(tokens, repeat(idx))
        return tokens

    hash_iter = chain(*map(hash_token, count(), sequences))
    hashes = np.fromiter(hash_iter, dtype=[('hash', np.int64), ('pos', np.int64)])

    hashes = np.sort(hashes, order='hash')
    grouped_hashes = np.split(hashes, np.unique(hashes['hash'], return_index=True)[1][1:])
    grouped_seqs = [[sequences[i] for i in c['pos']] for c in filter(lambda g: g.shape[0] > 1, grouped_hashes)]
    # make groups of unique sequences
    grouped_seqs = [list(set(g)) for g in grouped_seqs]
    # filter out singletons and return
    return list(filter(lambda g: len(g) > 1, grouped_seqs))

    # sorted_hashes = np.sort(hashes, order='hash')
    # grouped_hashes = np.split(sorted_hashes, np.unique(sorted_hashes['hash'], return_index=True)[1][1:])
    # return [[sequences[i] for i in c['seq']] for c in filter(lambda g: g.shape[0] > 1, grouped_hashes)]


def levenshtein_cluster(sequences):
    '''
    Take a list of sequences and cluster them by levenshtein distance 1. This is a very slow algorithm and should only be used for small datasets because it is O(n^2).
    For indels, the algorithm will only consider insertions and deletions of length 1. For substitutions, the algorithm will only consider substitutions of length 1.
    The clusters are returned as a list of lists of strings, where all strings in a list are within levenshtein distance 1 of each other.
    The clusters do not include insertions, only deletion of larger sequence length. This is because the algorithm is designed for following network graph construction, which generate all relevant edges of 1-edit type over all possible sequences combinations.
    :param sequences: the list of sequences to cluster (must be strings)
    :return: the clusters
    '''

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


def levenshtein_cluster2(sequences: List[str]) -> List[List[str]]:
    """
    Clusters a list of strings based on their Levenshtein distance of 1.

    This function takes a list of strings as its input and clusters them based on their Levenshtein distance of 1.
    The algorithm considers only 1-deletions for generating 1-edit sequences.

    The clusters are returned as a list of lists, where each list contains strings that are within a Levenshtein distance of 1 from each other.
    The clusters do not include insertions, only deletions of larger sequence length.

    Parameters:
        sequences (List[str]): The list of strings to cluster.

    Returns:
        List[List[str]]: The list of clusters.
    """
    # Generate all 1-deletion sequences
    deletion_sequences = starmap(combinations, zip(sequences, repeat(len(sequences))))
    # Hash the deletion sequences
    deletion_hashes = [hash(s) for s in chain(*deletion_sequences)]
    # Find the sequences that hash to the same value
    deletion_clusters = defaultdict(set)
    for i, h in enumerate(deletion_hashes):
        deletion_clusters[h].add(deletion_sequences[i])
    # Add the original sequences to the clusters
    for i, s in enumerate(sequences):
        deletion_clusters[hash(s)].add(s)
    # Filter out singletons
    # deletion_clusters = [list(c) for c in deletion_clusters.values() if len(c) > 1]
    # Filter out clusters that contain sequences of different lengths
    # deletion_clusters = [c for c in deletion_clusters if len(set(len(s) for s in c)) == 1]
    return deletion_clusters


# seq_gen = partial(sequence_generator, min_len=2, max_len=6)
# x = list(seq_gen(1000))
import random

# x = sequence_generator(100000, min_len=8, max_len=20)
x = rand_rep(1000000)
# print()
# for each sequence in mutated_sequences, randomly mutate by replacing a 1 random other character with a random character
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
mutated_sequences = random.sample(x, len(x) // 10)
for seq in mutated_sequences:
    i = random.randint(0, len(seq) - 1)
    aa = seq[i]
    replaced = set(AA)
    replaced.remove(aa)
    mutated = seq[:i] + random.choice(list(replaced)) + seq[i + 1:]
    x.append(mutated)
random.shuffle(x)

# runtime
import time
start = time.time()
clusters = levenshtein_cluster(x)
end = time.time()
print(end - start)

start = time.time()
clusters = levenshtein_cluster2(x)
end = time.time()
print(end - start)
