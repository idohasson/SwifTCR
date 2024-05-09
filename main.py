import csv
from collections import defaultdict
from itertools import combinations

import networkx as nx
import pandas as pd
import tqdm
from more_itertools import map_reduce

# Generate all (length-1) substrings of the current sequence and hash them
def generate_hashed_substrings(sequence):
    for substring in combinations(sequence, len(sequence) - 1):
        yield hash(substring), sequence


input_file_path = "./tcr.txt"  # Column name in the file to extract sequences from
column_name = 'aaSeqCDR3'  # Delimiter used in the input file (e.g., ',' for comma-separated)
default_delimiter = '\t'  # Output file to write the clustered sequences
cluster_file_path = 'dist_one_clusters.csv'  # Output file to write the edge list
edge_list_file_path = "dist_one_edge_list.csv"  # Header for the edge list CSV file
edge_list_header = ['Seq1', 'Seq2']


# # Read sequences from the input file
# with open(input_file_path) as f:
#     n_rows = sum(1 for line in f)
#
# sequences = []
# read_bar = tqdm.tqdm(total=n_rows, desc="Reading sequences")
# for chunk in pd.read_csv(input_file_path, sep=default_delimiter, usecols=[column_name], chunksize=512):
#     sequences.extend(chunk[column_name].to_list())
#     read_bar.update(len(chunk))
# read_bar.close()
#
# # Filter out invalid sequences (non-string or length <= 1)
# valid_sequences = filter(lambda s: isinstance(s, str) and len(s) > 1, sequences)
# # Group sequences by their lengths
# sequences = map_reduce(valid_sequences, keyfunc=len, reducefunc=lambda x: tuple(set(x)))

# FOR TESTING
import random  # Generate all (length-1) substrings of the current sequence
sequences = map_reduce((''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=l))
                        for l in range(2, 6) for _ in range(10 ** l)),
                       keyfunc=len, reducefunc=lambda x: tuple(set(x)))


# Write the header row for the edge list CSV file
with open(edge_list_file_path, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(edge_list_header)

write_edge_list_bar = tqdm.tqdm(total=sum(len(s) * l for l, s in sequences.items()), desc="Writing pairs")
for current_length in sorted(sequences.keys(), reverse=True):

    sequences_of_current_length = tuple(map(generate_hashed_substrings, sequences[current_length]))

    sequences_of_following_length = []
    if current_length - 1 in sequences.keys():
        sequences_of_following_length = tuple(map(lambda s: hash(tuple(s)), sequences[current_length - 1]))

    for i in range(current_length - 1, -1, -1):
        # Create a dictionary to store sequences grouped by their hashed subsequences
        sequences_by_hashed_substrings = defaultdict(list)

        for hashed_subsequence, sequence in map(next, sequences_of_current_length):
            sequences_by_hashed_substrings[hashed_subsequence].append(sequence)

        # If sequences of following length exist, add them to the dictionary as indel sequences
        if current_length - 1 in sequences.keys():
            for index, hashed_sequence in enumerate(sequences_of_following_length):
                if hashed_sequence in sequences_by_hashed_substrings:
                    sequences_by_hashed_substrings[hashed_sequence].append(sequences[current_length - 1][index])

        # Filter out groups with only one sequence
        sequences_by_hashed_substrings = [v for k, v in sequences_by_hashed_substrings.items() if len(v) > 1]

        edge_list = []
        for group in sequences_by_hashed_substrings:  # Sort sequences and generate edges
            if current_length - 1 in sequences.keys() and len(group[-1]) == current_length - 1:
                group[:-1] = sorted(group[:-1])
                edge_list.extend(combinations(group[:-1], 2))
                edge_list.extend([(s, group[-1]) for s in group[:-1] if i == 0 or (i > 0 and s[i] != s[i - 1])])
            else:
                group = sorted(group)
                edge_list.extend(combinations(group, 2))

        # Write the edge list to a CSV file
        with open(edge_list_file_path, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(edge_list)

        write_edge_list_bar.update(len(sequences[current_length]))
    sequences.pop(current_length)
write_edge_list_bar.close()

# Read the edge list from the CSV file
edges_list = []
with open(edge_list_file_path) as f:
    reader = csv.reader(f)
    next(reader)
    edges_list += [tuple(row) for row in reader]

# Create a new graph
graph = nx.Graph()
# Add edges from the edge list to the graph
graph.add_edges_from(edges_list)
# # Find connected components in the graph
# connected_components = list(nx.connected_components(graph))
# Find cliques in the graph
clique_list = ("|".join(clique) for clique in nx.find_cliques(graph))
# Write the edge list to a CSV file
with open(cluster_file_path, 'w', newline='') as file:
    for line in tqdm.tqdm(clique_list, desc="Writing clusters"):
        file.write(line + '\n')


# ###########################################        Test        #####################################################
# # !!!!!!!!!! remove the last line in the main loop sequences.pop(current_length) or re-read the sequences !!!!!!!!!!
# ####################################################################################################################
# from collections import Counter
# from Levenshtein import distance, editops
# from operator import eq
# gt = defaultdict(set)
# for sequence in set.union(*map(set, sequences.values())):
#     for ss in combinations(sequence, len(sequence) - 1):
#         gt[hash(ss)].add(sequence)
#     gt[hash(tuple(sequence))].add(sequence)
# gt = set(pair for v in filter(lambda v: len(v) > 1, gt.values())
#          for pair in combinations(sorted(v), 2)
#          if len(pair[0]) != len(pair[1]) or sum(1 for p in zip(*pair) if p[0] != p[1]) == 1)
# edge_list = set(map(tuple, map(sorted, edges_list)))
# print("Edges Test:", "Correct" if len(edges_list) == len(gt) and
#                                   len(edge_list) == len(gt) and
#                                   not gt - edge_list | edge_list - gt
#                                   and all(v == 1 for v in Counter(edges_list).values()) else "Incorrect")
#
# print("Clusters Test: ", "Correct" if all(distance(*pair) == 1 and eq(*editops(*pair)[0][1:])
#                                           for clique in clique_list for pair in
#                                           combinations(clique.split("|"), 2)) else "Incorrect")
