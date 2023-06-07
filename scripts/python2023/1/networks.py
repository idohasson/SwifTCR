import networkx as nx

from SwifTCR.linkage import get_edges


def get_sequence_network(sequences, method='lev'):
    """
    Returns a NetworkX graph object representing the network of valid edges between sequences
    using the specified distance metric.
    :param sequences: list of strings, the sequences to generate edges for
    :param method: string, the distance metric to use, one of 'hamming', 'lv', or 'damerau'
    :return: NetworkX graph object, the network of valid edges between sequences
    """
    # Get the valid edges between sequences using the specified distance metric
    edges = get_edges(sequences, method)

    # Return a NetworkX graph object representing the edges
    return nx.from_edgelist(edges)

