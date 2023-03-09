from itertools import *

from more_itertools import *

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
    return prepend(tuple(seq), combinations(seq, len(seq)-1))

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



def subsequence_generator(sequence):
    """
    Find subsequences with the length of the sequences minus one.
    """
    for i in range(len(sequence)):
        yield subsequence(sequence, i)




def subsequence_hash(sequence):
    for skip in range(len(sequence)):
        yield hash(tuple(aa for i, aa in enumerate(sequence) if i != skip)), skip


def hash_sequences(sequences):
    return np.vstack([np.array((ss_hash, skip_pos, seq_i))
                      for seq_i, seq in enumerate(sequences)
                      for ss_hash, skip_pos in subsequence_hash(seq)])