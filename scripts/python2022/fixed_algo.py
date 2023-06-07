import csv
import sys
import networkx as nx

def get_data(file_path, nrow, read_col=1, min_size=5):
    sequence_collection = dict()
    with open(file_path, "r") as csvfile:
        data_reader = csv.reader(csvfile)
        next(data_reader)  # the header row
        i = 0
        for row in data_reader:
            aa_seq = row[read_col]
            if isinstance(aa_seq, str):
                seq_len = len(aa_seq)
                if seq_len > min_size:
                    if seq_len not in sequence_collection:
                        sequence_collection[seq_len] = set()
                    sequence_collection[seq_len].add(aa_seq)
                    i += 1
                    if i == nrow:
                        break
    return sequence_collection


def sub_seq_table(i, sequence_collection, no_gaps=True):
    subs_table = dict()
    for seq in sequence_collection:
        sub_seq = seq[:i] + seq[i + 1:]
        subs_table.setdefault(sub_seq, set()).add(seq)

    return {k: v for k, v in subs_table.items() if
            len(v) > 1} if no_gaps else subs_table  # remove singletons to avoid gaps


def tables_by_len(sequence_l, sequence_collection):
    len_tables = dict()
    for index in range(sequence_l):
        len_tables[index] = sub_seq_table(index, sequence_collection)
    return len_tables


def get_hash_table(sequence_len_dict):
    hash_table = dict()
    for seq_len, sequence_set in sequence_len_dict.items():
        hash_table[seq_len] = tables_by_len(seq_len, sequence_set)
    return hash_table


# NO GAPS
def main(input_path, output_path, read_number):
    dataset = get_data(input_path, read_number)
    tables = get_hash_table(dataset)
    pass


if __name__ == '__main__':
    i_path = sys.argv[1]
    o_path = sys.argv[2]
    row_n = 100000
    # row_n = sys.argv[3]
    # main(i_path, o_path, int(row_n))
    main(i_path, o_path, row_n)
