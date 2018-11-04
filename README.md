# kcollections
[![Build Status](https://travis-ci.com/masakistan/kcollections.svg?token=oruFeF6Jkw9aGsjG6xUW&branch=master)](https://travis-ci.com/masakistan/kcollections)

## About
`kcollections` is a python-like `dict` and `set` that efficiently store kmers as keys or elements of the collection.
It is designed for building/prototyping Bioinformatics tools that rely on kmers but where the included `dict` and `set` consume too much memory for use.

It implements the [Bloom Filter Trie](https://github.com/GuillaumeHolley/BloomFilterTrie) algorithm.
This implementation differs from Guillaume et al. by allowing kmers of arbitrary size and by providing a generic dictionary/map data structure for associating arbitrary values with kmers.

`kcollections` is currently only available for Python version 2.7.

## Installation
The recommended way of installing `kcollections` is:

`pip install kcollections`

Alternatively, you can build from source.

### Build from source
`kcollections` must be built using `clang` due to a bug that exists in [`gcc`](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=36566).
Prerequisites include:

  - [jemalloc](http://jemalloc.net/)
  - [pybind11](https://github.com/pybind/pybind11)
  
These prerequisites are included or retrieved automatically using the `cmake` or `setup.py` build tools.

To build and install the python module from source:

```bash
git clone https://github.com/masakistan/kcollections.git
cd kcollections
pip install .
```

## Performance
`kcollections` is quite a bit slower than the `dict` or `set` but is much more memory-efficient.
We measured memory usage and running time using `/usr/bin/time -v` on a `Intel(R) Xeon(R) E5-2650v4 @2.20GHz` with 256 GB RAM.
27mers used for testing were taken from the human genome, no repetitive kmers appear in our dataset providing a worst case scenario where no insertions or queries are pruned before traversing the entire data structure.

### Memory Usage

|# kmers|`kset`|`set`|
|-------|------|-----|
|1 million|25.32 MB|105.82 MB|
|10 million|63.74 MB|906.96 MB|
|100 million|0.56 GB|11.98 GB|
|500 million|2.42 GB|48.54 GB|
|1 billion|4.43 GB|97.07 GB|
|1.5 billion|6.44 GB|191.61 GB|
|2 billion|8.44 GB||
|2.4 billion|10.08 GB|220.06 GB|

![Figure of memory usage](./memory_fig.png)

### Insertion Time
Times are `h:mm:ss`.

|# kmers|`kset`|`set`|
|-------|------|-----|
|1 million|0:00:03|0:00:01|
|10 million|0:00:36|0:00:13|
|100 million|0:07:45|0:02:20|
|500 million|0:41:13|0:13:10|
|1 billion|1:27:07|0:33:01|
|1.5 billion|2:11:11|0:57:42|
|2 billion|2:52:22|1:13:36|
|2.4 billion|3:24:58|1:30:26|

![Figure of insertion time](./insert_fig.png)

### Querying Time
Times are `h:mm:ss`.

|# kmers|`kset`|`set`|
|-------|------|-----|
|1 million|0:00:03|0:00:01|
|10 million|0:00:33|0:00:11|
|100 million|0:06:17|0:02:02|
|500 million|0:36:33|0:12:26|
|1 billion|1:18:40|0:30:24|
|1.5 billion|2:02:06|0:55:00|
|2 billion|2:44:30|1:15:59|
|2.4 billion|3:16:52|2:02:14|

![Figure of query time](./query_fig.png)

## Example Usage

### Using Kdict

```python
import kcollections
kd = kcollections.Kdict(27)

# insertion and value assignment
kd['AAACTGTCTTCCTTTATTTGTTCAGGG'] = 'banana phone'
assert kd['AAACTGTCTTCCTTTATTTGTTCAGGG'] == 'banana phone'

# iteration
for kmer, val in kd.iteritems():
    print kmer, val
```

### Using Kset

```python
import kcollections
ks = kcollections.Kset(27)

# insertion
ks.add('AAACTGTCTTCCTTTATTTGTTCAGGG')
assert 'AAACTGTCTTCCTTTATTTGTTCAGGG' in ks

# iteration
for kmer in ks:
    print kmer
```

### Read Mapper and Assembler
An example read mapping algorithm and assembler are provided using `kcollections` in the `applications` directory.

## Citation

## Acknowledgements
`kcollections` was built at the Computational Science Laboratory at Brigham Young University by Stanley Fujimoto (@masakistan) and Cole Lyman (@colelyman).

Funding was provided by the Utah NASA Space Grant Consortium and the BYU Graduate Research Fellowship.
