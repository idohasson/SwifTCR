def get_hashes(combinations_sets):
    merged_sets = itertools.chain(*combinations_sets)
    hash_sets = map(hash, merged_sets)
    return hash_sets


def read_file(file_path):
    # n = 0
    # pp = int(row_number / 100)
    with open(file_path, "r") as csvfile:
        datareader = csv.reader(csvfile)
        next(datareader)  # yield the header row
        for row in datareader:
            yield row[1]

            # n += 1
            # if not n % pp:
            #     clear_output(wait=True)
            #     print(np.round(100 * n / row_number, 2), "%", end="", flush=True)
            # if n == row_number:
            #     clear_output(wait=True)
            #     return


# def get_seq_combinations(full_sequences, seq_len):
#     return set(itertools.combinations(full_sequences, seq_len - 1))

def get_seq_combinations(full_sequences):
    return set(itertools.combinations(full_sequences, len(full_sequences) - 1))


def get_hash_comb(sequence_reader):
    sil = []
    hl = []
    n = 0
    for seq in sequence_reader:
        if seq is not None:
            comb = get_seq_combinations(seq)
            h_combs = list(map(hash, comb))
            hl += h_combs

            combs_n = len(h_combs)
            sil += [seq] * combs_n
            n += combs_n

    return hl, sil


def get_comb(sequence_reader):
    sil = []
    hl = []
    n = 0

    comb = defaultdict(list)
    full_seq = defaultdict(list)
    # y = defaultdict(list)
    for seq in sequence_reader:
        l = len(seq)
        # comb = get_seq_combinations(seq, l)
        comb[l].append(get_seq_combinations(seq, l))
        full_seq[l].append(seq)
        # y[l].append(len(comb))
        # seq_holder.append(get_seq_combinations(seq))

        # comb = get_seq_combinations(seq)
        # h_combs = list(map(hash, comb))
        # hl += h_combs
        #
        # combs_n = len(h_combs)
        # sil += [seq] * combs_n
        # n += combs_n

    return full_seq, comb


def sorted_indices(to_sort):
    return np.argsort(np.array(to_sort))


def get_cliques(hl, sorted_i):
    clique_s = []
    clique = set()
    for si1, si2 in itertools.zip_longest(sorted_i[:-1], sorted_i[1:]):
        if hl[si1] == hl[si2]:
            clique.add(si1)
            clique.add(si2)
        elif clique:
            clique_s.append(clique.copy())
            clique.clear()
    return clique_s


def get_edges(from_clique, seq_i):
    em = defaultdict(list)
    for cs in from_clique:
        si_cs = [seq_i[oi] for oi in cs]
        e = list(itertools.permutations(si_cs, 2))
        for u, v in e:
            em[u].append(v)
    return em


def find_clusters(edges):
    visited = dict.fromkeys(edges.keys(), False)
    queue = [list(edges.keys())[0]]
    w_cbuffer = [[]]
    while queue:
        v = queue.pop()
        if not visited[v]:
            visited[v] = True
            queue += [u for u in edges[v] if not visited[u]]
            w_cbuffer[-1].append(v)

        if not queue:
            if not all(visited.values()):
                for i, x in visited.items():
                    if not x:
                        queue.append(i)
                        break
                for j in w_cbuffer[-1]:
                    del visited[j]
                w_cbuffer.append([])
            # cid += 1
    return w_cbuffer


def write_clusters(w_path, w_clusters):
    with open(w_path, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(["seq_id", "cluster_id"])
        # write the data
        for cid, w_cluster in enumerate(w_clusters):
            for wc in w_cluster:
                writer.writerow([wc, cid])


def main(input_path, output_path):
    seq_reader = read_file(input_path)

    hash_list, seq_i_list = get_hash_comb(seq_reader)

    sorted_hashes = sorted_indices(hash_list)

    clique_sets = get_cliques(hash_list, sorted_hashes)

    if not clique_sets:
        print("No clusters")
        return

    edge_map = get_edges(clique_sets, seq_i_list)

    clusters = find_clusters(edge_map)

    write_clusters(output_path, clusters)


if __name__ == '__main__':
    import csv
    import sys
    import itertools
    import numpy as np
    from collections import defaultdict
    from IPython.display import clear_output

    i_path = sys.argv[1]
    o_path = sys.argv[2]
    main(i_path, o_path)
