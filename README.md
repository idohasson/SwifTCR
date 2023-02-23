# SwiftTCR: A Quick-and-Simple Way to Discover T Cell Repertoire Patterns
SwifTCR is a Python package provides a relatively quick way to find all sequences with shared single-only-replacement, only-deletion or both, creating clusters of CDR3s highly similar and construct a network graph. The tool, designed to enable researchers to quickly identify clusters in a relatively small data-set, is also available for large-scale data collection thanks its implementation of Spark, making it available for big-data immunological analysis.

## Installation

Install from [PyPI](https://pypi.org/project/SwifTCR/) by typing into your terminal or command-line:

``` bash
pip install swiftcr
```

## Examples

#### Cluster sequence lists

```python
from swifTCR import get_clusters 

seq_list = [['ABC', 'AXC', 'ABBC'], ['ABB', 'AC']]
get_clusters(seq_list, repalcement_only=False)

```

#### Generate random CDR3 and create a networkX graph:

```python
from SwifTCR import rand_rep, get_network

random_cdr3_list = rand_rep(seq_n=100, min_len=3, max_len=8)
# Single amino acid replacement network graph
G = get_network(random_cdr3_list)

```

#### Cluster from multiple MiCXR output files in a directory with spark

```python
from swifTCR import spark_cluster_file

spark_cluster('path/to/dir', 
                file_type='mixcr', 
                filename='.*TRCB.txt',
                out='path/to/save/spark/csv/')

```


## Documentation

[Documentation](https://linktodocumentation)


## License

[MIT](https://choosealicense.com/licenses/mit/)

