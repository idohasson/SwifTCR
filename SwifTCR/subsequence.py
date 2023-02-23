from itertools import combinations, count
from more_itertools import flatten, chunked, prepend

'''
Create a substring for each string in the input collection
including all possible versions of the missing letter. The
logic behind this strategy is that obtaining mutual same-
substing-combination amongst two or more sequences implies 
that they must be precisely one replacement apart, meaning 
that the hamming edit distance is one.
'''

def hashSubStr(sequence):
    return map(hash, combinations(sequence, len(sequence)-1))

def subsequence(seq):
    '''Return L-1 length subsequences of the input string.'''
    return combinations(seq, len(seq)-1)

def sub_and_sequence(seq):
    '''Return iterable of tuples L-1 length sub-strings and string's tuple.'''
    return prepend(tuple(seq), subsequence(seq))

def hash_subsequences(subs):
    '''
    Return an iterable containing the hash values of each sub-sequence of 
    the iterable sub-string.
    '''
    return zip(map(hash, subs), count())

def hash_key_gen(sequence_list, subsequences_only=False, chunk_size=None):
    '''
    Returns a hash value iterable of all sub-sequences of the given sequences.
    A sub-sequence is a sequence that does not contain one of the letters.
    '''
    ss_func = subsequence if subsequences_only else sub_and_sequence
    ss = map(ss_func, sequence_list)
    hash_ss = map(hash_subsequences, ss)
    return chunked(flatten(hash_ss), chunk_size)


from itertools import combinations, starmap
from operator import add
from more_itertools import flatten

def subsequence(seq):
    '''Return L-1 length subsequences of the input string.'''
    return combinations(seq, len(seq)-1)
def subseq_hash(seq):
    '''Return hashed sub-sequenceces of the input string, each pair with its index.'''
    return map(hash, subsequence(seq))
def indexed_hash(seq):
    '''return an iterable of every sub-sequence's hash value summed with its index of creation order'''
    return starmap(add, enumerate(subseq_hash(seq)))
def indexed_hash_iter(sequence_list):
    '''return list of long integers of hash values summed with the subsequence's index'''
    return flatten(map(indexed_hash, sequence_list))
