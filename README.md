# SwifTCR: A Quick-and-Simple Approach to Discover T-Cell Receptor Repertoire Patterns
SwifTCR is a Python package that identifies clusters of highly similar CDR3 sequences with a Hamming or Levenshtein edit distance of 1. It can be used to quickly identify clusters in datasets by detecting shared single-only-replacement, only-deletion, or both.

## Examples

#### Cluster sequence list by a single edit distance:

```python
from SwifTCR import find_clusters

find_clusters(['ABC', 'AXC', 'ABBC', 'ABB', 'AC'])
# Output:
# [{'ABC', 'AXC'}, {'ABB', 'ABC'}]

find_clusters(['ABC', 'AXC', 'ABBC', 'ABB', 'AC'], dist_type="levenshtein")
# Output:
# [(2, {'ABC', 'AC', 'AXC'}),
#  (1, {'ABB', 'ABBC'}),
#  (1, {'ABB', 'ABC'}),
#  (2, {'ABBC', 'ABC'}),
#  (3, {'ABBC', 'ABC'})]
```

## License

[MIT](https://choosealicense.com/licenses/mit/)

