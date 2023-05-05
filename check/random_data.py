import random
import scipy.stats as ss

__all__ = ['rand_rep']


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
    # generate random lengths
    lengths = normal_dist_int(min_len, max_len, seq_n, 2)
    # generate random sequences
    sequences = [''.join(random.choices(alphabet, k=length)) for length in lengths]
    return sequences

def normal_dist_int(l, h, n, sigma):
    """
    Return a list of length n of random integers from the normal distribution
    with mean and standard deviation as the mean of l and h and sigma,
    respectively. The values are clipped to be between l and h.
    """
    mean = (l + h) / 2
    std = (h - l) / (2 * sigma)
    return ss.truncnorm((l - mean) / std, (h - mean) / std, loc=mean, scale=std).rvs(n).astype(int)


AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


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