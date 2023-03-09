import operator
from itertools import starmap

from more_itertools import pairwise

def is_hamming_distance_one(seq1, seq2):
    """
    Returns True if the Hamming distance between two sequences is 1, False otherwise.
    :param seq1: first sequence
    :param seq2: second sequence
    :return: boolean indicating if the Hamming distance between the sequences is 1
    """
    return sum(starmap(operator.ne, zip(seq1, seq2))) == 1

def is_levenstein_distance_one(seq1, seq2):
    """
    Returns True if the Levenstein distance between two sequences is 1, False otherwise.
    :param seq1: first sequence
    :param seq2: second sequence
    :return: boolean indicating if the Levenstein distance between the sequences is 1
    """

    if abs(len(seq1) - len(seq2)) == 1:
        return True

    return is_hamming_distance_one(seq1, seq2)


def is_damerau_distance_one(seq1, seq2):
    """
    Returns True if the Damerau-Levenstein distance between two sequences is 1, False otherwise.
    :param seq1: first sequence
    :param seq2: second sequence
    :return: boolean indicating if the Damerau-Levenstein distance between the sequences is 1
    """
    if is_levenstein_distance_one(seq1, seq2):
        return True

    return sum((aa11 != aa21 and aa12 != aa22) and (aa11 == aa22 and aa12 == aa21) for
               (aa11, aa12), (aa21, aa22) in pairwise(zip(seq1, seq2))) == 1


#
# def is_valid_pairing(substitution=True, deletion=True, transposition=True):
#     """
#     Returns a function that checks if two sequences are valid pairings, based on the selected pairing methods.
#     :param substitution: boolean, whether to use substitution as a pairing method
#     :param deletion: boolean, whether to use deletion as a pairing method
#     :param transposition: boolean, whether to use transposition as a pairing method
#     :return: function that takes two sequences as input and returns a boolean indicating if they are valid pairings
#     """
#     # Define filter functions for each pairing method
#     filters = {
#         'substitution': lambda x: sum(starmap(operator.ne, zip(x[0], x[1]))) == 1,
#         'deletion': lambda x: abs(len(x[0]) - len(x[1])) == 1,
#         'transposition': lambda x: sum((aa11 != aa21 and aa12 != aa22) and (aa11 == aa22 and aa12 == aa21) for
#                                           (aa11, aa12), (aa21, aa22) in pairwise(zip(x[0], x[1]))) == 1
#     }
#
#     # Raise an error if no valid pairing method was selected
#     if not any([substitution, deletion, transposition]):
#         raise ValueError('At least one pairing method must be selected. Please set one or more of the arguments '
#                          '"substitution", "deletion", or "transposition" to True.')
#
#     selected_filters = []
#     if substitution:
#         selected_filters.append(filters['substitution'])
#     if deletion:
#         selected_filters.append(filters['deletion'])
#     if transposition:
#         selected_filters.append(filters['transposition'])
#
#         # Return a function that checks if two sequences are valid pairings
#
#     return lambda x: any(f(x) for f in selected_filters)



