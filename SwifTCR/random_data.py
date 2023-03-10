import random
from itertools import filterfalse

import scipy.stats as ss
from scipy.stats import norm

__all__ = ['rand_rep']

from more_itertools import random_combination

import random
from math import exp, sqrt, pi

import random
import math

def normal_dist_choices(values, k, sigma=1):
    n = len(values)
    mu = (n - 1) / 2  # Mean index of the values list
    weights = [math.exp(-(i - mu)**2 / (2*sigma**2)) for i in range(n)]
    counts = [math.ceil(w * k) for w in weights]
    indices = random.sample(range(n), k=k, counts=counts)
    return [values[i] for i in indices]

values = [1, 2, 3, 4, 5]
indices = normal_dist_choices(values, k=3, sigma=1)
print(indices)  # e.g. [2, 4, 1]



values = [1, 2, 3, 4, 5]
result = normal_dist_choices(values, 3)
print(result)



def normal_dist_int(l, h, n, sigma=2):
    """
    Return a list of length n of random integers from the normal distribution
    with mean and standard deviation as the mean of l and h and sigma,
    respectively. The values are clipped to be between l and h.
    """
    mean = (l + h) / 2
    std = (h - l) / (2 * sigma)
    return ss.truncnorm((l - mean) / std, (h - mean) / std, loc=mean, scale=std).rvs(n).astype(int)

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
    # lengths = ss.truncnorm((min_len - mean) / std, (max_len - mean) / std, loc=mean, scale=std).rvs(seq_n).astype(int)
    # lengths = normal_dist_int(min_len, max_len, seq_n)
    # generate random sequences
    # lengths = random.choices(list(range(min_len, max_len + 1)), weights=list(norm()) , k=seq_n)
    sequences = [''.join(random.choices(alphabet, k=length)) for length in lengths]
    # Ensure that the data have pairwise 1 edit distance of all kinds
    # and some more with more than 1 edit distance.
    mutation = {
        'hamming': lambda seq, pos, aa: seq[:pos] + aa + seq[pos + 1:],
        'levenstein': lambda seq, pos, aa: random_combination(seq, len(seq)-1),
        'levenstein-del': lambda seq, pos, aa: random_combination(seq, len(seq)-1),
        'levenstein-insert': lambda seq, pos, aa: seq[:pos] + aa + seq[pos:],
        'damerau': lambda seq, pos, aa: seq[:pos] + seq[pos+1] + seq[pos] + seq[pos + 2:],
    }

    mutated_sequences = []
    for seq in random.choices(list(filter(lambda x: len(x) > 1, sequences)), k=len(sequences)//10):
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





# AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def generate_aa(l):
    """
    Generate random amino acids sequence.
    """
    aa = "".join(random.choice(AA, size=l))
    return aa


def sequence_generator(n, min_len=6, max_len=17):
    """
    Generate random N amino acids sequences with random length.
    """
    lengths = random.choice(range(min_len, max_len + 1), size=n)
    sequences = [generate_aa(l) for l in lengths]
    return sequences

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