import operator
from itertools import starmap

from more_itertools import pairwise


def is_valid_pairing(substitution=True, deletion=True, transposition=True):
    """
    Returns a function that checks if two sequences are valid pairings, based on the selected pairing methods.
    :param substitution: boolean, whether to use substitution as a pairing method
    :param deletion: boolean, whether to use deletion as a pairing method
    :param transposition: boolean, whether to use transposition as a pairing method
    :return: function that takes two sequences as input and returns a boolean indicating if they are valid pairings
    """
    # Define filter functions for each pairing method
    filters = {
        'substitution': lambda a, b: sum(starmap(operator.ne, zip(a, b))) == 1,
        'deletion': lambda a, b: abs(len(a) - len(b)) == 1,
        'transposition': lambda a, b: sum((aa11 != aa21 and aa12 != aa22) and (aa11 == aa22 and aa12 == aa21) for
                                          (aa11, aa12), (aa21, aa22) in pairwise(zip(a, b))) == 1
    }

    # Select filter functions based on the input arguments
    selected_filters = [filters[key] for key in filters.keys() if locals()[key]]

    # Raise an error if no valid pairing method was selected
    if not selected_filters:
        raise ValueError('At least one pairing method must be selected. Please set one or more of the arguments '
                         '"substitution", "deletion", or "transposition" to True.')

    # Return a function that checks if two sequences are valid pairings based on the selected filter functions
    return lambda a, b: any(filter_func(a, b) for filter_func in selected_filters)

