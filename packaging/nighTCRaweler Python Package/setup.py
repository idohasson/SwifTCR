from setuptools import setup, find_packages

VERSION = '0.4.15'
DESCRIPTION = 'Cluster and construct netwotk graph of highly similar sequences'
LONG_DESCRIPTION = 'SwifTCR is a Python package provides a relatively quick way to find all sequences with shared single-only-replacement, only-deletion or both, creating clusters of CDR3s highly similar and construct a network graph. The tool, designed to enable researchers to quickly identify clusters in a relatively small data-set, is also available for large-scale data collection thanks its implementation of Spark, making it available for big-data immunological analysis.'
REQUIREMENTS = ["more_itertools", "scipy", "numpy", "pandas", "networkx"],

setup(
        name="SwifTCR", 
        version=VERSION,
        author="ido",
        author_email="idoh@systemsbiomed.org",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url='https://github.com/systemsbiomed/SwifTCR',
        packages=find_packages(),
        install_requires=REQUIREMENTS,
        keywords=['cluster', 'TCR', 'edit', 'distance'],
)