from itertools import filterfalse, compress

import numpy as np
import scipy.stats as ss
from scipy.stats import norm

__all__ = ['rand_rep']

from more_itertools import random_combination
from functools import partial
import random
import math

AA = list('ACDEFGHIKLMNPQRSTVWY')
NT = list('ACGT')

def rand_rep(seq_n, min_len=8, max_len=18, char_n=20, alphabet='ACDEFGHIKLMNPQRSTVWY'):
    '''
    Generate list of CDR3 amino acid sequences

    Generate strings with distribiutions of the repertoire's size, sequences' 
    size and amino acid types py position that are more or less similar to 
    real data samples of TCR repertoire.


    Parameters
    ----------
    seq_n : int
        Number of unique sequences
    min_len : int, optional
        Minimum sequence size to be generated, by default 8
    max_len : int, optional
        Maximum sequence size to be generated, by default 18

    Returns
    -------
    list
        Strings of amino acid sequences
    '''
    if isinstance(char_n, int):
        alphabet = alphabet[:char_n]

    sequence_lengths = normal_dist_int(min_len, max_len, seq_n, 2)
    return ["".join(random.choices(alphabet, k=seq_size)) for seq_size in sequence_lengths]


def normal_dist_int(l, h, n, sigma):
    lengths = np.arange(l, h + 1)
    d = lengths.size / 2
    x = np.arange(-d, d, 1)
    xU, xL = x + 1, x
    prob = ss.norm.cdf(xU, scale=sigma) - ss.norm.cdf(xL, scale=sigma)
    prob = prob / prob.sum()
    nums = np.random.choice(lengths, size=n, p=prob)
    return nums


def normal_dist_choices(values, k, sigma=1):
    n = len(values)
    mu = (n - 1) / 2  # Mean index of the values list
    weights = [math.exp(-(i - mu) ** 2 / (2 * sigma ** 2)) for i in range(n)]
    counts = [math.ceil(w * k) for w in weights]
    indices = random.sample(range(n), k=k, counts=counts)
    return [values[i] for i in indices]




def rand_rep(seq_n, min_len=8, max_len=18, char_n=20, alphabet='ACDEFGHIKLMNPQRSTVWY'):
    '''
    Generate list of CDR3 amino acid sequences

    Generate strings with distribiutions of the repertoire's size, sequences' 
    size and amino acid types py position that are more or less similar to 
    real data samples of TCR repertoire.


    Parameters
    ----------
    seq_n : int
        Number of unique sequences
    min_len : int, optional
        Minimum sequence size to be generated, by default 8
    max_len : int, optional
        Maximum sequence size to be generated, by default 18

    Returns
    -------
    list
        Strings of amino acid sequences
    '''
    # sigma = 2
    # # generate random lengths
    # mean = (min_len + max_len) / 2
    # std = (max_len - min_len) / (2 * sigma)
    # # generate random sequences
    # lengths = ss.truncnorm((min_len - mean) / std, (max_len - mean) / std, loc=mean, scale=std).rvs(seq_n).astype(int)
    # lengths = normal_dist_int(min_len, max_len, seq_n)
    lengths = random.choices(list(range(min_len, max_len + 1)), weights=list(norm()), k=seq_n)
    sequences = [''.join(random.choices(alphabet, k=length)) for length in lengths]
    # Ensure that the data have pairwise 1 edit distance of all kinds
    # and some more with more than 1 edit distance.
    mutation = {
        'hamming': lambda seq, pos, aa: seq[:pos] + aa + seq[pos + 1:],
        'levenstein': lambda seq, pos, aa: random_combination(seq, len(seq) - 1),
        'levenstein-del': lambda seq, pos, aa: random_combination(seq, len(seq) - 1),
        'levenstein-insert': lambda seq, pos, aa: seq[:pos] + aa + seq[pos:],
        'damerau': lambda seq, pos, aa: seq[:pos] + seq[pos + 1] + seq[pos] + seq[pos + 2:],
    }

    mutated_sequences = []
    for seq in random.choices(list(filter(lambda x: len(x) > 1, sequences)), k=len(sequences) // 10):
        mutation_type = random.choice(['hamming', 'levenstein-del', 'levenstein-insert', 'damerau'])
        if mutation_type == 'damerau':
            pos = random.randint(0, len(seq) - 2)
        elif mutation_type == 'levenstein-del':
            pos = random.randint(0, len(seq) - 1)
        else:
            pos = random.randint(0, len(seq) - 1)

        if mutation_type == 'hamming':

            aa = random.choice(list(filterfalse(lambda x: x == seq[pos], alphabet)))
        elif mutation_type == 'levenstein' or mutation_type == 'levenstein-insert':
            aa = random.choice(alphabet)
        else:
            aa = ''

        mut_func = mutation[mutation_type]
        mutated = mut_func(seq, pos, aa)
        mutated_sequences.append(mutated)

    sequences.extend(mutated_sequences)
    random.shuffle(sequences)
    return sequences


def random_aa_seq(seq_size):
    """
    Generate random amino acids sequence.
    """

    return "".join(random.choices(population=('A', 'C', 'D', 'E', 'F',
                                              'G', 'H', 'I', 'K', 'L',
                                              'M', 'N', 'P', 'Q', 'R',
                                              'S', 'T', 'V', 'W', 'Y'),
                                  k=seq_size))


import random
from itertools import chain


def apply_mutation(seq, mutation_type, alphabet='ACDEFGHIKLMNPQRSTVWY'):
    """
    Apply a random mutation of a given type to a sequence and return the original and mutated sequences.
    """
    i = random.randint(0, len(seq) - 1)
    original_seq = seq

    yield original_seq

    if mutation_type == 'substitution':
        mutation = random.choice(list(set(alphabet) - {seq[i]}))
        mutated_seq = ''.join(chain(seq[:i], mutation, seq[i + 1:]))
    elif mutation_type == 'deletion':
        mutated_seq = ''.join(chain(seq[:i], seq[i + 1:]))
    elif mutation_type == 'insertion':
        mutation = random.choice(alphabet)
        mutated_seq = ''.join(chain(seq[:i], mutation, seq[i:]))
    elif mutation_type == 'transposition':
        j = random.randint(0, len(seq) - 1)
        mutated_seq = ''.join(chain(seq[:i], seq[j], seq[i + 1:j], seq[i], seq[j + 1:]))
    else:
        raise ValueError('Invalid mutation type.')

    yield mutated_seq



def sequence_generator(seq_num, min_len=1, max_len=5):
    """
    Generate random N amino acids sequences with random length.
    """
    # generate random lengths at normal distribution range of min_len and max_len (inclusive)
    lengths = random.choices(range(min_len, max_len + 1), k=seq_num)
    # generate random sequences with random lengths
    sequences = map(random_aa_seq, lengths)

    for seq in sequences:
        yield seq

# import random
# x = sequence_generator(1000, min_len=10, max_len=15)
# # print()
# # for each sequence in mutated_sequences, randomly mutate by replacing a 1 random other character with a random character
# AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
# mutated_sequences = random.sample(x, len(x) // 10)
# for seq in mutated_sequences:
#     i = random.randint(0, len(seq) - 1)
#     aa = seq[i]
#     replaced = set(AA)
#     replaced.remove(aa)
#     mutated = seq[:i] + random.choice(list(replaced)) + seq[i + 1:]
#     x.append(mutated)
# random.shuffle(x)
