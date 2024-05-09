import csv
import tqdm
from collections import defaultdict
from itertools import combinations
from more_itertools import map_reduce
import pandas as pd


# Generate all (length-1) substrings of the current sequence and hash them
def generate_hashed_substrings(sequence):
    for substring in combinations(sequence, len(sequence) - 1):
        yield hash(substring), sequence


input_file_path = "./tcr.csv"  # Column name in the file to extract sequences from
column_name = 'aaSeqCDR3'  # Delimiter used in the input file (e.g., ',' for comma-separated)
default_delimiter = '\t'  # Output file to write the clustered sequences
output_file_path = 'clustered_dist_one.txt'  # Output file to write the edge list
edge_list_file_path = "edge_list_dist_one.csv"  # Header for the edge list CSV file
edge_list_header = ['Seq1', 'Seq2']

# Write the header row for the edge list CSV file
with open(edge_list_file_path, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(edge_list_header)

# Read sequences from the input file
with open(input_file_path) as f:
    n_rows = sum(1 for line in f)

sequences = []
read_bar = tqdm.tqdm(total=n_rows, desc="Reading sequences")
for chunk in pd.read_csv(input_file_path, sep=default_delimiter, usecols=[column_name], chunksize=512):
    sequences.extend(chunk[column_name].to_list())
    read_bar.update(len(chunk))
read_bar.close()

# Filter out invalid sequences (non-string or length <= 1)
valid_sequences = filter(lambda s: isinstance(s, str) and len(s) > 1, sequences)
# Group sequences by their lengths
sequences = map_reduce(valid_sequences, keyfunc=len, reducefunc=set)

# # Generate all (length-1) substrings of the current sequence
# import random
# sequences = (''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=l)) for l in range(2,6) for _ in range(2 ** l))
# sequences = map_reduce(sequences, len, reducefunc=set)


write_bar = tqdm.tqdm(total=sum(len(s) * l for l, s in sequences.items()), desc="Writing pairs")
for current_length in sorted(sequences.keys(), reverse=True):

    sequences_of_current_length = list(map(generate_hashed_substrings, sequences[current_length]))

    sequences_of_following_length = []
    if current_length - 1 in sequences.keys():
        sequences_of_following_length = list(map(lambda s: (hash(tuple(s)), s), sequences[current_length - 1]))

    for i in range(current_length - 1, -1, -1):

        # Create a dictionary to store sequences grouped by their hashed subsequences
        sequences_by_hashed_substrings = defaultdict(list)

        for hashed_subsequence, sequence in map(next, sequences_of_current_length):
            sequences_by_hashed_substrings[hashed_subsequence].append(sequence)

        # If sequences of following length exist, add them to the dictionary as indel sequences
        if current_length - 1 in sequences.keys():
            for hashed_sequence, sequence in sequences_of_following_length:
                if hashed_sequence in sequences_by_hashed_substrings:
                    sequences_by_hashed_substrings[hashed_sequence].append(sequence)

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

        write_bar.update(len(sequences[current_length]))
write_bar.close()

# # Test
# gt = defaultdict(set)
# for sequence in set.union(*sequences.values()):
#     for ss in combinations(sequence, len(sequence) - 1):
#         gt[hash(ss)].add(sequence)
#     gt[hash(tuple(sequence))].add(sequence)
# gt = set(pair for v in filter(lambda v: len(v) > 1, gt.values())
#          for pair in combinations(sorted(v), 2)
#          if len(pair[0]) != len(pair[1]) or sum(1 for p in zip(*pair) if p[0] != p[1]) == 1)
#
# # Read the edge list from the CSV file
# edges_from_file = []
# with open(edge_list_file_path) as f:
#     reader = csv.reader(f)
#     next(reader)
#     edges_from_file += [tuple(row) for row in reader]
#
# test_edge_list = set(map(tuple, map(sorted, edges_from_file)))
# print("Test:", "Correct" if len(edges_from_file) == len(gt) and
#                             len(test_edge_list) == len(gt) and
#                             not gt - test_edge_list | test_edge_list - gt
#                             and all(v == 1 for v in Counter(edges_from_file).values()) else "Incorrect")
