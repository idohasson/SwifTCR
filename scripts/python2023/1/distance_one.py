from SwifTCR.clustering import get_clusters
from SwifTCR.linkage import get_edges


def get_distance_one_clusters(sequences, method='lev'):
    '''
    Returns a list of clusters of sequences that are distance 1 apart from each other.
    The clustering is done using the zero/one deletion approach, which groups sequences
    sharing a common subsequence with a single deletion.
    :param sequences: list of strings, the sequences to cluster
    :param method: string, the clustering method to use, one of 'hamming', 'lev', or 'damerau'
    :return: list of lists of strings, the clusters of sequences that are distance 1 apart from each other
    '''
    # Get the clusters of sequences using the specified clustering method
    clusters = get_clusters(sequences, method)

    # Return the list of clusters
    return clusters


def get_distance_one_edges(sequences, method='lev'):
    '''
    Generates a generator of valid edges for a set of sequences that are distance 1 apart from each other.
    The edges are generated using the zero/one deletion approach, which clusters sequences sharing a common
    subsequence with a single deletion. An edit distance validation method suited for the specified edit-distance
    method is then applied to each sequence pair.
    :param sequences: list of strings, the sequences to generate edges for
    :param method: string, the distance metric to use, one of 'hamming', 'lev', or 'damerau'
    :return: generator of tuples of strings, the valid edges between sequences that are distance 1 apart from each other
    '''
    # Get the clusters of sequences that are distance 1 apart from each other using the specified clustering method
    clusters = get_clusters(sequences, method)

    # Generate the valid edges between sequences in the clusters using the specified distance metric
    edges = get_edges(clusters, method)

    # Return the generator of valid edges
    return edges

