import string
from collections import defaultdict
from functools import partial
from itertools import combinations, chain, groupby, product
from operator import itemgetter

import numpy as np


def get_tokens(seq, depth=1):
    return {seq: [''.join(s) for s in powerset(seq) if len(seq) - len(s) <= depth]}

def join_dicts(dicts):
    return {k: v for d in dicts for k, v in d.items()}

def get_corpus(seq, depth=1):
    token_func = partial(get_tokens, depth=depth)
    return join_dicts(map(token_func, seq))

# function that make e
def get_tokens_dict(tokens):
    d=defaultdict(set)
    for s, t in tokens.items():
        for ss in t:
            d[ss].add(s)
        d[s].add(s)
    return d

# frumtion that make e
def get_pairs(d):
    return set(chain(*map(lambda g: combinations(sorted(g), 2), d.values())))

from more_itertools import bucket, powerset


def read_words(filename):
    with open(filename) as f:
        return {word.strip().lower() for word in f}


def create_graph(words):
    graph = {word: set() for word in words}
    for word in words:
        for i, char in product(range(len(word) + 1), string.ascii_lowercase):
            new_word = word[:i] + char + word[i:]
            if new_word in graph and new_word != word:
                graph[word].add(new_word)
                graph[new_word].add(word)
        for i, char in product(range(len(word)), string.ascii_lowercase):
            new_word = word[:i] + char + word[i+1:]
            if new_word in graph and new_word != word:
                graph[word].add(new_word)
                graph[new_word].add(word)
    return graph


def find_groups(graph):
    groups = []
    seen = set()
    for word in graph:
        if word not in seen:
            group = set(bucket(powerset(graph[word]), 2))
            seen.update(group)
            groups.append(group)
    return groups


def find_largest_group(groups):
    return max(groups, key=len)


def main():
    words = read_words('/usr/share/dict/words')
    graph = create_graph(words)
    groups = find_groups(graph)
    largest_group = find_largest_group(groups)
    print('Largest group:', largest_group)



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



def levenshtein_cluster2(sequences):
    """
    Given a list of sequences, cluster them by the Levenshtein distance 1.
    @param sequences - the list of sequences to cluster
    @return The list of clusters
    """
    hashes = np.array([(hash(s[:i] + s[i + 1:]), i, id)
                       for id, s in enumerate(sequences)
                       for i in range(len(s) + 1)])
    # sort by the hash and index keeping the order of the id
    sorted_hashes = np.sort(hashes, order=['f0', 'f1'])
    # get the unique hashes and the index of the first occurence
    unique_hash, index = np.unique(sorted_hashes[:, 0], return_index=True)
    # get the index of the first occurrence of the next unique hash and subtract the index of the first occurrence of the current unique hash
    # this gives the number of times the hash occurs in the array
    # if the difference between the current and previous unique hash is greater than 1
    # then the current hash is a group of sequences that are not within 1 Levenshtein distance of each other
    # so we split the array at this point and move on to the next group of sequences
    index = [index[i] for i in range(1, len(unique_hash)) if abs(unique_hash[i] - unique_hash[i - 1]) > 1]
    # split the array at the indices
    clusters = np.split(sorted_hashes, index, axis=0)
    # for each cluster, get the unique sequences and the index of the first occurrence of each sequence
    # the index of the first occurrence of each sequence is the index of the sequence in the original array

    pass
    return





# test run time of levenshtein_cluster


# import timeit
# print(timeit.timeit(lambda: len(levenshtein_cluster(x)), number=1))
# print(timeit.timeit(lambda: len(levenshtein_cluster2(x)), number=1))

# d=get_tokens_dict(get_corpus(x))
# print(get_pairs(d))


CHARACTERS = list("abcdefghijklmnopqrstuvwxyz")
graph = {}

# ~ 240 000 words
for word in x:
    graph[word] = set()

for val in graph:
    node = graph[val]
    for char in CHARACTERS:
        for i in range(len(val) + 1):
            word = val[:i] + char + val[i:]
            if word in graph and word != val:
                node2 = graph[word]
                node.add(node2)
                node2.add(node)
        for i in range(len(val)):
            word = val[:i] + char + val[i+1:]
            if word in graph and word != val:
                node2 = graph[word]
                node.add(node2)
                node2.add(node)
