import operator
from itertools import starmap, tee, islice
from more_itertools import pairwise, exactly_n, quantify, locate, zip_offset


# def n_substitutions(seq1, seq2, n=1):
#     """
#     Returns the number of substitutions between two sequences.
#     :param seq1: first sequence
#     :param seq2: second sequence
#     :return: number of substitutions between the sequences
#     """
#
#     if len(seq1) != len(seq2):
#         raise ValueError('Sequences must be of the same length.')
#
#     return exactly_n(zip(seq1, seq2), n=n, predicate=lambda x: x[0] != x[1])
#
#
# def n_deletions(seq1, seq2, n=1):
#     """
#     Returns the number of deletions between two sequences.
#     :param seq1: first sequence
#     :param seq2: second sequence
#     :return: number of deletions between the sequences
#     """
#
#     return abs(len(seq1) - len(seq2)) == n
#
# # def n_transpositions(seq1, seq2, n=1):
# #
# #     if len(seq1) != len(seq2):
# #         raise ValueError('Sequences must be of the same length.')
# #
# #     cross, diff = tee(zip(seq1, seq2))
# #
# #     adjacent_crossing = starmap(lambda x, y: operator.eq(x[1], y[0]) and operator.eq(x[0], y[1]), pairwise(cross))
# #
# #     cross_match = map(operator.and_, starmap(operator.ne, diff), adjacent_crossing)
# #
# #     return quantify(cross_match) <= n
#
# def n_transpositions(seq1, seq2, n=1):
#
#     if len(seq1) != len(seq2):
#         raise ValueError('Sequences must be of the same length.')
#
#     if n_substitutions(seq1, seq2, n+1) != 2:
#         return False
#
#     loc = locate()
#
#
#     unmatch_pairs = locate(zip(seq1, seq2), lambda x: operator.ne(*x))
#
#     adjacent_pairs = ((aa11 != aa21 and aa12 != aa22) and (aa11 == aa22 and aa12 == aa21) for
#                (aa11, aa12), (aa21, aa22) in pairwise(zip(seq1, seq2)))
#
#     return unmatch_pairs


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

    diff = map(operator.ne, seq1, seq2)
    consecutive_diff = starmap(operator.and_, pairwise(diff))
    one_consecutive_diff = exactly_n(consecutive_diff, n=1)
    return one_consecutive_diff

    # if exactly_n(n_diff, n=2, predicate=lambda x: operator.ne(*x)):
    #     return False
    #
    # return sum((aa11 != aa21 and aa12 != aa22) and (aa11 == aa22 and aa12 == aa21) for
    #            (aa11, aa12), (aa21, aa22) in pairwise(cross)) == 1

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
