from itertools import combinations, chain

# from SwifTCR.pairing import is_valid_pairing


def get_skip_links(clusters, valid_func):
    """
    Returns a generator of valid edges for a set of clusters using a specified pairing method.
    :param clusters: iterable of iterables of strings, the clusters of sequences to pair
    :param valid_func: function, the pairing method to use
    :return: generator of tuples of strings, the valid edges between sequences in the clusters
    """
    # Define a lambda function to get all edges between sequences in a cluster
    cluster_edges = lambda cluster: combinations(sorted(cluster), 2)

    # Define a lambda function to check if a pair of sequences is a valid pairing using the specified function
    valid_pair = lambda pair: valid_func(pair[0], pair[1])

    # Define a lambda function to filter the valid edges for a cluster
    filter_func = lambda cluster: filter(valid_pair, cluster_edges(cluster))

    # Get the valid edges for each cluster and combine them into a single iterator
    valid_edges = map(filter_func, clusters)
    edges_found = chain.from_iterable(valid_edges)

    # Return the iterator of valid edges
    return edges_found


def get_edges(clusters, method='lev'):
    """
    Returns a generator of valid edges for a set of clusters using the specified distance metric.
    :param clusters: iterable of iterables of strings, the clusters of sequences to pair
    :param method: string, the distance metric to use, one of 'hamming', 'lev', or 'damerau'
    :return: generator of tuples of strings, the valid edges between sequences in the clusters
    """
    # Ensure that the specified distance metric is valid
    assert method in ['hamming', 'lev', 'damerau']

    # Define the distance metric function for each method
    dist_method = {
        'hamming': is_valid_pairing(substitution=True, deletion=False, transposition=False),
        'lev': is_valid_pairing(substitution=True, deletion=True, transposition=False),
        'damerau': is_valid_pairing(substitution=True, deletion=True, transposition=True)
    }

    # Select the appropriate distance metric function based on the specified method
    valid_func = dist_method[method]

    # Get the valid edges for the clusters using the selected distance metric function
    edges = get_skip_links(clusters, valid_func)

    # Return the iterator of valid edges
    return edges
