from collections import defaultdict
from itertools import combinations
from networkx import Graph, connected_components

def find_connected_groups(sequences, k=8):
    # Extract all unique k-mers from the sequences
    kmers = set()
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmers.add(seq[i:i+k])

    # Construct a graph where each node represents a k-mer, and each edge
    # represents a mutation (i.e., substitution, insertion, or deletion) between two k-mers
    graph = Graph()
    for kmer in kmers:
        graph.add_node(kmer)
    for kmer1, kmer2 in combinations(kmers, 2):
        if calculate_hamming_distance(kmer1, kmer2) == 1:
            graph.add_edge(kmer1, kmer2)

    # Find connected components in the graph, which represent sets of sequences that are reachable
    # from one another through a sequence of mutations
    groups = []
    for component in connected_components(graph):
        group = set()
        for kmer in component:
            for seq in sequences:
                if kmer in seq:
                    group.add(seq)
        groups.append(group)

    return groups

def calculate_hamming_distance(s1, s2):
    # Calculate the Hamming distance between two strings
    return sum(1 for a, b in zip(s1, s2) if a != b)

# Example usage:
sequences = ["ACGTACGT", "ACGTTAGT", "ACGTAGCT", "TGTAGCTT"]
groups = find_connected_groups(sequences)
print(groups)
