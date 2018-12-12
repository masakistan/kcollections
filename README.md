# kcollections
[![Build Status](https://travis-ci.com/masakistan/kcollections.svg?token=oruFeF6Jkw9aGsjG6xUW&branch=master)](https://travis-ci.com/masakistan/kcollections)

## About
`kcollections` is a python-like `dict` and `set` that efficiently store kmers as keys or elements of the collection.
It is designed for building/prototyping Bioinformatics tools that rely on kmers but where the included `dict` and `set` consume too much memory for use.

It implements the [Bloom Filter Trie](https://github.com/GuillaumeHolley/BloomFilterTrie) algorithm.
This implementation differs from Guillaume et al. by allowing kmers of arbitrary size and by providing a generic dictionary/map data structure for associating arbitrary values with kmers.

`kcollections` is currently only available for Python version 2.7.

## Installation
We provide some pre-compiled binaries for python 2/3 and Linux and MacOS:

```bash
pip install kcollections
```

Alternatively, you can build from source.

### Build from source
Prerequisites include:

  - [jemalloc](http://jemalloc.net/)
  - [pybind11](https://github.com/pybind/pybind11)

These prerequisites are included or retrieved automatically using the `cmake` or `setup.py` build tools.

To build and install the python module from source:

```bash
git clone https://github.com/masakistan/kcollections.git
cd kcollections

# python 3
python3 setup.py bdist_wheel
pip3 install dist/*.whl

# python 2
python2 setup.py bdist_wheel
pip2 install dist/*.whl
```

## Performance
`kcollections` is quite a bit slower than the `dict` or `set` but is much more memory-efficient.
We measured memory usage and running time using `/usr/bin/time -v` on a`Intel(R)
Xeon(R) E5-2650v4 @2.20GHz` with 256 GB RAM using `Clang`.
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
|2 billion|8.44 GB|194.14 GB|
|2.4 billion|10.08 GB|220.06 GB|

![Figure of memory usage](./memory_fig.png)

### Insertion Time
Insertion time comparisons using built-in Python set, `kcollections` serial and parallel insert.

![Figure of insertion time](./insert_time.png)


## Example Usage

### Using Kdict

```python
import kcollections
kd = kcollections.Kdict(27)

# insertion and value assignment
kd['AAACTGTCTTCCTTTATTTGTTCAGGG'] = 'banana'
kd['AAACTGTCTTCCTTTATTTGTTCAGGT'] = 'phone'
assert kd['AAACTGTCTTCCTTTATTTGTTCAGGG'] == 'banana'
assert kd['AAACTGTCTTCCTTTATTTGTTCAGGT'] == 'phone'

# iteration
for kmer, val in kd.iteritems():
    print kmer, val

# removal
del kd['AAACTGTCTTCCTTTATTTGTTCAGGT']
```

### Using Kset

Kmers can be added one at a time with `add`, but the fastest way to add kmers to a set is
to add an DNA sequence using `add_seq`.
Faster still, use `parallel_add_seq` for multithreaded inserts.

```python
import kcollections
ks = kcollections.Kset(27)

# add single kmer
ks.add('AAACTGTCTTCCTTTATTTGTTCAGGG')

# sequence insertion
seq = 'AAACTGTCTTCCTTTATTTGTTCAGGGATCGTGTCAGTA'
ks.add_seq(seq, len(seq))

assert 'AAACTGTCTTCCTTTATTTGTTCAGGG' in ks

# multithreaded sequence insertion
# nthreads must be a power of 2.
# nthreads = 4 or 16 work well
ks.parallel_add_init(16)
ks.parallel_add_seq(seq, len(seq))
ks.parallel_add_join()

# iteration
for kmer in ks:
    print kmer
```

### Using Kcounter

`Kcounter` is an implementation of the Python collection's
[Counter](https://docs.python.org/2/library/collections.html#collections.Counter),
but the keys must be kmers, of course!
Like `Kdict`, kmers can be added to `Kcounter` one at a time, but the
fastest ways to add kmers to a set is to add an DNA sequence using `add_seq` (or
`parallel_add_seq` for multithreaded inserts).

``` python
from kcollections import Kcounter
kc = Kcounter(27)

# add single kmer
kc['AAACTGTCTTCCTTTATTTGTTCAGGG'] += 1

# sequence insertion
seq = 'AAACTGTCTTCCTTTATTTGTTCAGGGATCGTGTCAGTA'
kc.add_seq(seq, len(seq))

assert 'AAACTGTCTTCCTTTATTTGTTCAGGG' in kc

# multithreaded sequence insertion
# nthreads must be a power of 2.
# nthreads = 4 or 16 work well
kc.parallel_add_init(4)
kc.parallel_add_seq(seq, len(seq))
kc.parallel_add_join()

# iteration
for kmer, count in kc.iteritems():
    print kmer, count
```

### Read Mapper and Assembler
An example read mapping algorithm and assembler are provided using `kcollections` in the `applications` directory.

## Citation

## Acknowledgements
`kcollections` was built at the Computational Science Laboratory at Brigham Young University by Stanley Fujimoto (@masakistan) and Cole Lyman (@colelyman).

Funding was provided by the Utah NASA Space Grant Consortium and the BYU Graduate Research Fellowship.
