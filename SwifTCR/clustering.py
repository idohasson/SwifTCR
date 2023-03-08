from collections import defaultdict
from functools import partial
from itertools import combinations


def get_skip_groups(sequences, zero_deletion=True):
    """
    Returns a list of clusters of sequences with a single deletion or no deletion.
    :param sequences: list of strings, the sequences to cluster
    :param zero_deletion: boolean, whether to include clusters of sequences with no deletion
    :return: list of lists of strings, the clusters of sequences
    """
    # Create a defaultdict to store the clusters
    cluster_dict = defaultdict(list)

    # For each sequence:
    for seq in sequences:
        # Create a subsequence by removing one residue at a time,
        for subseq in combinations(seq, len(seq) - 1):
            # Add the original sequence to the cluster corresponding to that subsequence
            cluster_dict[''.join(subseq)].append(seq)
        # If zero_deletion is True, also add the original sequence to its own cluster
        if zero_deletion:
            cluster_dict[seq].append(seq)

    # Return the list of clusters
    return cluster_dict.values()


def get_clusters(sequences, method='lev'):
    """
    Returns a list of clusters of similar sequences using the specified method.
    :param sequences: list of strings, the sequences to cluster
    :param method: string, the clustering method to use, one of 'hamming', 'lev', or 'damerau'
    :return: list of lists of strings, the clusters of sequences
    """
    # Ensure that the specified method is valid
    assert method in ['hamming', 'lev', 'damerau']

    # Define the cluster functions for each method
    cluster_method = {
        'hamming': partial(get_skip_groups, zero_deletion=False),
        'lev': partial(get_skip_groups, zero_deletion=True),
        'damerau': partial(get_skip_groups, zero_deletion=True)
    }

    # Select the appropriate cluster function based on the specified method
    cluster_func = cluster_method[method]

    # Cluster the sequences using the selected function
    clusters = cluster_func(sequences)

    # Return the list of clusters
    return clusters

