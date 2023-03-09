from itertools import combinations, chain, count, repeat, filterfalse
from more_itertools import all_equal, pairwise, prepend, split_when
from collections import defaultdict
from operator import itemgetter
from functools import partial
import numpy as np

__all__ = ['find_clusters']


# def find_matching_sequences(sequences):
#     subseq_to_seqs = defaultdict(set)
#     for seq in sequences:
#         for l in range(1, len(seq)):
#             for subseq in combinations(seq, l):
#                 subseq_to_seqs[''.join(sorted(subseq))].add(seq)
#
#     matching_seqs = [list(seq_set) for seq_set in subseq_to_seqs.values() if len(seq_set) > 1]
#     return matching_seqs


def get_hashes(sequence_list, prepend_sequence=False):

    if prepend_sequence:
        hash_sub_fn = partial(hash_subsequence, prepend_sequence=True)
    else:
        hash_sub_fn = hash_subsequence

    hashes = chain(*map(hash_sub_fn, sequence_list))

    return hashes



def group_sequences(sequence_list, prepend_sequence=False):

    hash_iterable = get_hashes(sequence_list, prepend_sequence=prepend_sequence)

    hash_arr = np.fromiter(hash_iterable, dtype=[('hash', 'i8'), ('pos', 'i8')])

    sorted_indices = np.lexsort((hash_arr['pos'], hash_arr['hash']))

    hash_groups = split_when(sorted_indices, lambda s1, s2: s1[0] != s2[0])

    return hash_groups


# def group_hashes(hashes):
#
#     hash_groups = split_when(hashes, lambda s1, s2: s1[0] != s2[0])
#
#     return hash_groups


# def filter_groups(hash_groups):
#
#     hash_groups = filterfalse(lambda g: len(g) == 1, hash_groups)
#
#     return hash_groups


def find_matching_sequences(sequences):
    subseq_to_seqs = defaultdict(set)
    for seq in sequences:
        for l in range(1, len(seq)):
            for subseq in combinations(seq, l):
                subseq_to_seqs[''.join(sorted(subseq))].add(seq)

    matching_seqs = [list(seq_set) for seq_set in subseq_to_seqs.values() if len(seq_set) > 1]
    return matching_seqs


def hash_subsequence(sequence, prepend_sequence=False):
    sub = combinations(sequence, len(sequence) - 1)

    if prepend_sequence:

        sub = prepend(tuple(sequence), sub)

        return zip(count(-1), map(hash, sub))

    else:

        return zip(count(), map(hash, sub))


# def get_hashes(sequence_list, prepend_sequence=False):
#     if prepend_sequence:
#         hash_sub_fn = partial(hash_subsequence, prepend_sequence=True)
#     else:
#         hash_sub_fn = hash_subsequence
#     hashes = chain(*map(hash_sub_fn, sequence_list))
#     return hashes


def sort_hash(hash_iterable):

    hash_arr = np.fromiter(hash_iterable, dtype=[('hash', 'i8'), ('pos', 'i8')])

    sorted_indices = np.lexsort((hash_arr['pos'], hash_arr['hash']))

    return sorted_indices


def group_sequences(sequence_list, prepend_sequence=False):

    hashes = get_hashes(sequence_list, prepend_sequence=prepend_sequence)

    sorted_indices = sort_hash(hashes)

    hash_groups = split_when(sorted_indices, lambda s1, s2: s1[0] != s2[0])

    return hash_groups


# def group_hashes(hashes):
#
#     hash_groups = split_when(hashes, lambda s1, s2: s1[0] != s2[0])
#
#     return hash_groups


# def filter_groups(hash_groups):
#
#     hash_groups = filterfalse(lambda g: len(g) == 1, hash_groups)
#
#     return hash_groups


def hamming(sequence_list):
    '''
    Given a list of sequences, return a list of sets of sequences that are all identical except for one character.

    Args:
    sequence_list: a list of sequences to be compared
    Returns:
    A list of sets of sequences that are all identical except for one character.
    '''

    def hash_subsequences(seq):
        hash_sub = map(hash, combinations(seq, len(seq) - 1))
        return zip(hash_sub, count(), repeat(seq))

    sub_keys = chain(*map(hash_subsequences, sequence_list))
    sorted_keys = sorted(sub_keys, key=lambda k: k[0] + k[1])

    hash_groups = split_when(sorted_keys, lambda s1, s2: s1[:2] != s2[:2])
    hash_groups = filterfalse(lambda g: len(g) == 1 or all_equal(map(itemgetter(2), g)), hash_groups)
    hash_groups = map(lambda g: set(map(itemgetter(2), g)), hash_groups)

    return hash_groups


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


def demarlio_levenshtein(sequence_list):
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

    def transposition(in_hash_group):
        return [((p1, p2), {s1, s2}) for (p1, c1), (p2, c2) in pairwise(in_hash_group)
                if p2 - p1 == 1 for s1 in c1 for s2 in c2 if s1[-p1] == s2[-p2] and s1[-p1] < s1[-p2]]

    def hash_clusters(groups):

        clusters = [(c[0][1], set(map(itemgetter(-1), c)))
                    for c in split_when(groups, lambda s1, s2: s1[1] != s2[1])]
        if clusters[0][0] == 0:
            tp_clusters = transposition(clusters[1:])
            clusters = [(p, c.union(clusters[0][1])) for p, c in clusters[1:]]
        else:
            tp_clusters = transposition(clusters)
        return filterfalse(lambda g: len(g[1]) == 1, chain(clusters, tp_clusters))

    sub_keys = chain(*map(hash_subsequences, sequence_list))
    sorted_keys = sorted(sub_keys, key=lambda k: k[0] + k[1])

    hash_groups = split_when(sorted_keys, lambda s1, s2: s1[0] != s2[0])
    hash_groups = filterfalse(lambda g: len(g) == 1, hash_groups)
    hash_groups = chain(*map(hash_clusters, hash_groups))
    return hash_groups


def find_clusters(sequence_list, dist_type="hamming", edge_list=False):
    '''
    The function find_clusters takes a list of sequences and a distance function.
    It returns a list of clusters, where each cluster is a list of sequences.

    Args:
    sequence_list: a list of sequences
    dist_type: The distance function to use.
    edge_list: If True, the function will return a list of edges. If False, it will return a list of nodes.
    Returns:
    A list of lists. Each list is a cluster of sequences.
    '''

    if dist_type == "hamming":
        dist_func = hamming
    elif dist_type == "levenshtein":
        dist_func = levenshtein
    elif dist_type == "demarlio-levenshtein":
        dist_func = demarlio_levenshtein
    else:
        return

    clusters = dist_func(sequence_list)
    return list(clusters)
