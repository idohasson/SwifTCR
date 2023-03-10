# Generate data incluse 1 edit distances of:
# - Hamming
# - Levenstein
# - Damerau-Levenstein
#

import random

from SwifTCR.random_data import rand_rep
def generate_data(n, l_min, l_max, alphabet):
    """
    Generate a list of random strings of the specified length and alphabet.
    :param n: number of strings to generate
    :param length: length of the strings to generate
    :param alphabet: alphabet to generate the strings from
    :return: list of random strings
    """
    random.choices(alphabet, k=n)
    data = rand_rep(n, l_min, l_max, alphabet)
    # Ensure that the data have pairwise 1 edit distance of all kinds
    # and some more with more than 1 edit distance.
    sample = random.sample(data, k=2)
    return sample