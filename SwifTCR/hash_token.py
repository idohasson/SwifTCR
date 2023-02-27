import random
import time
from collections import defaultdict
from itertools import chain, filterfalse
from itertools import combinations, count, repeat
from operator import itemgetter
from more_itertools import all_equal, pairwise, prepend, split_when


# create wrapper for timing functions
def timeit(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f'{func.__name__} took {end - start:.3f} seconds')
        return result

    return wrapper


def apply_combinations_dict(d, seq):
    return chain.from_iterable(d.get(y, set()) for y in combinations(seq, len(seq) - 1))


@timeit
def create_combinations_dict(sequences):
    d = defaultdict(set)
    for x, seq in enumerate(sequences):
        for y, z in zip(combinations(seq, len(seq) - 1), zip(count(), repeat(x))):
            d[y].add(z)
    return d


@timeit
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


@timeit
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


# test runtimes
sequences = [''.join(random.choices('ACGT', k=10)) for _ in range(100000)]
d = create_combinations_dict(sequences)
# ---------------------------- #
c1 = hamming(sequences)
# ---------------------------- #
c2 = levenshtein(sequences)
# ---------------------------- #
c3 = demarlio_levenshtein(sequences)
