
from SwifTCR.clustering import get_clusters
from SwifTCR.linkage import get_edges

import pandas as pd
def get_cluster_dataframe(sequences, method='lev'):
    """
    Returns a Pandas DataFrame representing the clusters of sequences using the specified clustering method.
    :param sequences: list of strings, the sequences to cluster
    :param method: string, the clustering method to use, one of 'hamming', 'lev', or 'damerau'
    :return: Pandas DataFrame, the clusters of sequences
    """
    # Get the clusters of sequences using the specified clustering method
    clusters = list(get_clusters(sequences, method))

    # Convert the clusters to a Pandas DataFrame
    cluster_df = pd.DataFrame({'cluster': clusters})

    # Split the cluster column into two separate columns for the sequences
    cluster_df[['seq1', 'seq2']] = pd.DataFrame(cluster_df['cluster'].tolist(), index=cluster_df.index)

    # Drop the original cluster column
    cluster_df = cluster_df.drop('cluster', axis=1)

    # Return the Pandas DataFrame
    return cluster_df


def get_cluster_dict(sequences, method='lev'):
    """
    Returns a dictionary of clusters of sequences using the specified clustering method,
    where each key is an integer index for the cluster and the value is a list of the sequences
    in the cluster.
    :param sequences: list of strings, the sequences to cluster
    :param method: string, the clustering method to use, one of 'hamming', 'lev', or 'damerau'
    :return: dictionary, the clusters of sequences
    """
    # Get the clusters of sequences using the specified clustering method
    clusters = list(get_clusters(sequences, method))

    # Convert the clusters to a dictionary
    cluster_dict = {i: cluster for i, cluster in enumerate(clusters)}

    # Return the dictionary of clusters
    return cluster_dict


def get_edge_dataframe(sequences, method='lev'):
    """
    Returns a Pandas DataFrame representing the valid edges between sequences that are distance 1 apart
    from each other using the specified distance metric.
    :param sequences: list of strings, the sequences to generate edges for
    :param method: string, the distance metric to use, one of 'hamming', 'lev', or 'damerau'
    :return: Pandas DataFrame, the valid edges between sequences that are distance 1 apart from each other
    """
    # Get the valid edges between sequences that are distance 1 apart from each other using the specified distance metric
    edges = get_edges(sequences, method)

    # Convert the edges to a Pandas DataFrame
    edge_df = pd.DataFrame(edges, columns=['seq1', 'seq2'])

    # Return the Pandas DataFrame
    return edge_df

