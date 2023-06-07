# AIRR-seq community format out headers
AIRR = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'junction_aa', 'duplicate_count',
            'rev_comp', 'productive', 'sequence_alignment', 'germline_alignment',
            'junction', 'v_cigar', 'd_cigar', 'j_cigar']

def prepare(sequence_list, min_len, max_len):
    '''
    Filter sequences by length.

    Given a list of sequences, return a list of sequences with length between min_len and max_len.
    '''
    if not isinstance(min_len, int) or min_len < 2:
        min_len = 2

    if not isinstance(max_len, int):
        return filter(lambda s: isinstance(s, str) and len(s) >= min_len, sequence_list)
    elif max_len <= min_len:
        return filter(lambda s: isinstance(s, str) and len(s) == min_len, sequence_list)
    else:
        return filter(lambda s: isinstance(s, str) and len(s) >= min_len and len(s) <= max_len, sequence_list)


def filter_size(sequences, min_len, max_len):
    if min_len < 2:
        min_len = 2
    if max_len is None:
        return [seq for seq in sequences if len(seq) >= min_len]
    return [seq for seq in sequences if min_len <= len(seq) <= max_len]

def prepare_cdr3(cdr3_lists, min_len=2, max_len=None):
    # filter each list by length
    cdr3_lists = [filter_size(seq, min_len, max_len) for seq in cdr3_lists]
    # concatenate all lists
    return [item for sublist in cdr3_lists for item in sublist]


# ---------------------------- I/O ----------------------------

def read_set_from_file(filename):
    """
    Extract a de-duped collection (set) of text from a file.
    Expected file format is one item per line.
    """
    collection = set()
    with open(filename, "r", encoding="utf-8") as file_:
        for line in file_:
            collection.add(line.rstrip())
    return collection


def get_file_extension(path, dot=True, lower=True):
    ext = os.path.splitext(path)[1]
    ext = ext if dot else ext[1:]
    return ext.lower() if lower else ext

def search_path(path):
    # check if path exists
    if not os.path.exists(path):
        raise FileNotFoundError(f'Path {path} does not exist')
    # check if path is a file
    if os.path.isfile(path):
        return [path]
    # check if path is a directory
    if os.path.isdir(path):
        # all files in directory
        files = os.listdir(path)
        # filter files by extension
        files = [os.path.join(path, file) for file in files if file.endswith('.csv') or
                 file.endswith('.tsv') or
                 file.endswith('.txt')]
        return files

def read_mixcr_cdr3(input_path, get_column='aaSeqCDR3', delim='\t'):
    '''
    Read mixcr output.

    Read and return list of the specified column data Path to the mixcr output from the given file / list of file / directory path.

    '''
    if isinstance(input_path, str):
        if os.path.isdir(input_path):
            input_path = [os.path.join(input_path, file) for file in os.listdir(input_path)]
        else:
            input_path = [input_path]

    def read_sequences(mixcr_file):
        if os.path.exists(mixcr_file) and os.path.isfile(mixcr_file):
            df = pd.read_csv(mixcr_file, header=0, sep=delim)
            if get_column in df.columns:
                return df[get_column].tolist()
        return []
    
    return [read_sequences(file_path) for file_path in input_path] if isinstance(input_path, list) else []




# ---------

import os




def search_path(path):
    """
     this function takes a path and returns a list of files in that path.
    @param path - the path to search for files
    @returns a list of files in the path
    """

    # check if path exists
    if not os.path.exists(path):
        raise FileNotFoundError(f'Path {path} does not exist')
    # check if path is a file
    if os.path.isfile(path):
        return [path]
    # check if path is a directory
    if os.path.isdir(path):
        # all files in directory
        files = os.listdir(path)
        # filter files by extension
        files = [os.path.join(path, file) for file in files if file.endswith('.csv') or
                 file.endswith('.tsv') or
                 file.endswith('.txt')]
        return files


def read_cdr3(file_paths, from_col=None, sep='\t'):
    """
    Read in the cdr3 sequences from the given file paths. If from_col is None, assume that the
    file is a tab-delimited file with the cdr3 sequences in the first column. If from_col is an int,
    assume that the file is a tab-delimited file with the cdr3 sequences in the column with that index.
    @param file_paths - the file paths to read from
    @param from_col - the column to read from
    @param sep - the separator to use
    @returns the cdr3 sequences
    """

    if from_col is None:
        from_col = 0

    if isinstance(from_col, str):
        read_func = lambda p: pd.read_csv(p, sep=sep, usecols=[from_col], header=0, dtype=str)
    elif isinstance(from_col, int):
        read_func = lambda p: pd.read_csv(p, sep=sep, index_col=from_col, header=0, dtype=str)
    else:
        raise ValueError(f'from_col must be str or int, not {type(from_col)}')
    # read all files in path_list
    return [read_func(file).aaSeqCDR3.to_list() for file in search_path(file_paths)]


def hash_search(cdr3_list, search_gaps=True, edge_list=True, out="pandas"):
    cdr3_list = filter(lambda aa: aa is not None and len(aa) >= 3, cdr3_list)
    clusters_found = distance_one(cdr3_list, search_gaps, edge_list)
    if clusters_found:
        clusters_found = [list(aa) + [i] for i, c in enumerate(clusters_found) for aa in c]
    else:
        return

    if out == "list":
        return clusters_found

    if out == "pandas" or out == "file":
        if edge_list:
            col_names = ['CDR3_1', 'CDR3_2', 'edit_position', 'cluster_id']
        else:
            col_names = ['CDR3', 'cluster_id']

    df_cluster = pd.DataFrame(clusters_found, columns=col_names)

    if out == "pandas":
        return df_cluster

    if out == "file":
        df_cluster.to_csv('output.tsv', sep='\t', index=False)



def read_pandas(file_path, field="aaSeqCDR3", sep="\t"):
    return pd.read_csv(file_path,
                       index_col=None,
                       header=0,
                       usecols=[field],
                       sep=sep,
                       dtype=str,
                       na_filter=False)



def cdr3_read(file_path, min_len=8, max_len=20, field="aaSeqCDR3", sep=","):
    df = pd.read_csv(file_path, index_col=None, header=0, usecols=[field], sep=sep, dtype=str,
                     na_filter=False).drop_duplicates()
    return {CDR3_length: group['aaSeqCDR3'].tolist() for CDR3_length, group in df.groupby(df.aaSeqCDR3.str.len()) if
            min_len <= CDR3_length <= max_len}



def cdr3_read(path):
    df = read_pandas(path)
    df = df.drop_duplicates()
    # TODO: filter invalid
    # by_len = df.groupby(df.aaSeqCDR3.str.len())
    # for idx, aa in df.groupby(df.aaSeqCDR3.str.len()):
    #     return idx, aa['aaSeqCDR3'].tolist()
        # return by_len.get_group(shorter_aa), by_len.get_group(aa)
    # frame['length'] = frame.aaSeqCDR3.str.len()
    # by_len = dict()
    # TODO: also the first len
    CDR3_by_length = dict()
    for CDR3_length, g in df.groupby(df.aaSeqCDR3.str.len()):
        CDR3_by_length[CDR3_length] = g['aaSeqCDR3'].tolist()
    return CDR3_by_length
            # yield CDR3_length, g['aaSeqCDR3'].tolist()

    # yield by_len
    # li = []
    # for filename in all_files:
    #     df = read_pandas(filename)
    #     li.append(df)


    # for a, b in combinations(range(len(li)), 2):
    #     frame = pd.concat([li[a], li[b]], axis=0, ignore_index=True).drop_duplicates()
    #     frame['length'] = frame.aaSeqCDR3.str.len()
    #
    #     by_len = dict()
    #     for name, group in frame.groupby('length'):
    #         by_len[name] = group['aaSeqCDR3'].tolist()
    #     yield by_len



def cdr3_preprocess(cdr3_list):
    """
    :param cdr3_list:
    :return:
    """
    by_length = {}
    for seq in cdr3_list:
        by_length.setdefault(len(seq), []).append(seq)
    for l_k, seqs in by_length.items():
        by_length[l_k] = list(dict.fromkeys(seqs))
    return by_length



def read_cdr3(file_path, field=None, sep="\t"):
    with open(file_path) as file:
        header = file.readline()
        col_index = header.strip().split(sep).index(field)
        # return list(line.strip().split(sep)[col_index].strip() for line in file)
        for line in file:
            yield line.strip().split(sep)[col_index].strip()



def read_pandas(file_path, field="aaSeqCDR3", sep=","):
    return pd.read_csv(file_path, index_col=None, header=0, usecols=[field],sep=sep, dtype=str, na_filter=False)

def cdr3_read(path, min_len=8, max_len=20):
    df = read_pandas(path).drop_duplicates()
    return {CDR3_length: g['aaSeqCDR3'].tolist() for CDR3_length, g in df.groupby(df.aaSeqCDR3.str.len()) if min_len <= CDR3_length <= max_len}

def read_chunk(file_path, n_seqs=None, chunk_size=None):
    for chunk in pd.read_csv(file_path, nrows=n_seqs, chunksize=chunk_size, header=0, usecols=["aaSeqCDR3"], squeeze=True):
        # filter sequences
        str_seq_list = filterfalse(lambda x: not isinstance(x, str), chunk)
        min_seq_len = filterfalse(lambda y: len(y) < 5, str_seq_list)
        yield min_seq_len



def cdr3_read(file_path, min_len=8, max_len=20, field="aaSeqCDR3", sep=","):
    df = pd.read_csv(file_path, index_col=None, header=0, usecols=[field], sep=sep, dtype=str,
                     na_filter=False).drop_duplicates()
    return {CDR3_length: group['aaSeqCDR3'].tolist() for CDR3_length, group in df.groupby(df.aaSeqCDR3.str.len()) if
            min_len <= CDR3_length <= max_len}


def read_cdr3(cdr3_path):
    df = pd.read_csv(cdr3_path, sep="\t", index_col=None, header=0, usecols=["aaSeqCDR3"])
    df = df.drop_duplicates()
    df = df.groupby(df.aaSeqCDR3.str.len())
    return df

