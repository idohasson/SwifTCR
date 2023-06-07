import itertools
import operator
import time
from functools import partial
import matplotlib.pyplot as plt
import numpy as np

from SwifTCR import rand_rep
from SwifTCR.cluster import get_skip_links
from SwifTCR.all.cluster import find_clusters


def runtime(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        # print(f'Function {func.__name__} took {end - start} seconds')
        return end - start
    return wrapper


hamm_func = partial(get_skip_links, method='hamming')
lev_func = partial(get_skip_links, method='lv')
dam_func = partial(get_skip_links, method='damerau')


# hamm_func2 = partial(dist_one_pairing2, method='hamming')
# lev_func2 = partial(dist_one_pairing2, method='lv')
# dam_func2 = partial(dist_one_pairing2, method='damerau')
# func_set2 = [hamm_func2, lev_func2, dam_func2]

find_hamming_clusters = partial(find_clusters, dist_type='hamming')
find_lev_clusters = partial(find_clusters, dist_type='levenshtein')
find_dam_clusters = partial(find_clusters, dist_type='demarlio-levenshtein')


times = {'Hamming': [], 'Levenshtein': [], 'Damerau-Levenshtein': []}

# sample_sizes = list(range(500_000, 1_000_001, 500_000))
sample_sizes = list(range(20000, 200001, 20000))
# sample_sizes = list(range(1000, 10000, 1000))
precentage = list(map(operator.truediv, itertools.accumulate(sample_sizes), itertools.repeat(sum(sample_sizes))))
func_set1 = [hamm_func, lev_func, dam_func]


# for n, prec in zip(sample_sizes, precentage):
#     sequences = rand_rep(n)
#     for dist_method, func1, func2 in zip(times.keys(), map(runtime, func_set1), map(runtime, func_set2)):
#         delta = func1(sequences) - func2(sequences)
#         times[dist_method].append(delta)
#     print(f'\rNumber of sequences: {n}, ({round(prec*100,1)}%)', end='')

for n, prec in zip(sample_sizes, precentage):
    sequences = rand_rep(n)
    for dist_method, func in zip(times.keys(), map(runtime, func_set1)):
        times[dist_method].append(func(sequences))
    print(f'\rNumber of sequences: {n}, ({round(prec*100,1)}%)', end='')

#     sequences = rand_rep(n)
#     for name, func in zip(['Hamming', 'Levenshtein', 'Damerau-Levenshtein'], map(runtime, [hamm_func, lev_func, dam_func])):
#         times[name].append(func(sequences))
#
#     print(f'\rNumber of sequences: {n}, ({round(prec*100,1)}%)', end='')


plt.style.use('seaborn-whitegrid')
plt.rcParams['figure.figsize'] = [10, 5]
plt.rcParams['figure.dpi'] = 100

for dist_type, t in times.items():
    plt.plot(sample_sizes, t, label=dist_type)

plt.legend()
plt.xlabel('Number of sequences')
plt.ylabel('Time (seconds)')
plt.title('Runtime of distance functions')

plt.show()

