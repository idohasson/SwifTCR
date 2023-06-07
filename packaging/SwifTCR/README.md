# SwifTCR: A Quick-and-Simple Approach to Discover T-Cell Receptor Repertoire Patterns
SwifTCR is a Python package that identifies clusters of highly similar CDR3 sequences with a Hamming or Levenshtein edit distance of 1. It can be used to quickly identify clusters in datasets by detecting shared single-only-replacement, only-deletion, or both.

## Examples

#### Cluster sequence list by a single edit distance:

```python
from SwifTCR import get_clusters

get_clusters(['ABC', 'AXC', 'ABBC', 'ABB', 'AC'])
# Output:
# [{'ABC', 'AXC'}, {'ABB', 'ABC'}]

get_clusters(['ABC', 'AXC', 'ABBC', 'ABB', 'AC'], method='lev')
# Output:
# [{'ABC', 'AXC', 'AC'}, {'ABB', 'ABBC'}]
```

#### Generate network of sequences connected by valid edges:

```python
from SwifTCR import get_sequence_network

sequences = ['ABC', 'AXC', 'ABBC', 'ABB', 'AC']
network = get_sequence_network(sequences, method='lev')

print(network.nodes())
# Output: ['ABC', 'AXC', 'ABBC', 'ABB', 'AC']

print(network.edges())
# Output: [('ABC', 'AXC'), ('ABC', 'AC'), ('AXC', 'AC'), ('ABB', 'ABBC'), ('ABB', 'ABC'), ('ABBC', 'ABC')]
```

#### Generate Pandas DataFrame of valid edges:

```python
from SwifTCR import get_edge_dataframe

sequences = ['ABC', 'AXC', 'ABBC', 'ABB', 'AC']
edge_df = get_edge_dataframe(sequences)

print(edge_df)
# Output:
#   seq1 seq2
# 0  ABC  AXC
# 1  ABC   AC
# 2  AXC   AC
# 3  ABB ABBC
# 4  ABB  ABC
# 5 ABBC  ABC
```

## License

[MIT](https://choosealicense.com/licenses/mit/)

