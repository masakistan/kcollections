# Kdict Merging for Parallel Operations

There are two differences between `Kset` and `Kdict` when performing parallel operations.
First, when initializing the data structure for parallel insertion `Kdict` requires a merging function be provided.
This tells the data structure what to do when there are multiple values being added for a single key.
Second, when using the `Kdict.parallel_add_seq` method a third argument is required that is an iterable.

## The Merging Function
The merging function must have two input parameters and one output.
For example:

```python
def merge_func(prev_val, new_val):
    return new_val
```

will overwrite the previous value associated with the kmer with the new value associated with the kmer.
The function can be much more complex depending on what you want to accomplish.

## Examples

### Kmer Counting
Though we provide `Kcounter`, kmer counting could be accomplished using `Kdict`.

``` python
dna = 'ACTGGTACTG'
kd = kcollections.Kdict(4)
kd.parallel_add_init(4, lambda prev_val, new_val: prev_val + new_val)
kd.parallel_add_seq(dna, len(dna), [1 for _ in range(len(dna))])
kd.parallel_add_join()

for kmer, val in kd.iteritems():
    print kmer, val
```

Output:
``` bash
GGTA 1
GTAC 1
CTGG 1
ACTG 2
TACT 1
TGGT 1
```

Here, the lambda function `lambda prev_val, new_val: prev_val + new_val` performs the counting.

### All Indices for Kmers
The following example stores the index of each kmer in a DNA string.
If a kmer occurs in multiple locations the value is converted into a list with all indices.

``` python
dna = 'ACTGGTACTG'

def merge_func(prev_val, new_val):
    if isinstance(prev_val, list):
        prev_obj.append(new_val)
        return prev_val
    else:
        return [prev_val, new_val]

kd = kcollections.Kdict(4)
kd.parallel_add_init(4, merge_func)
kd.parallel_add_seq(dna, len(dna), [i for i in range(len(dna))])
kd.parallel_add_join()

for kmer, val in kd.iteritems():
    print kmer, val
```

Output:

``` bash
GGTA 3
GTAC 4
CTGG 1
ACTG [0, 6]
TACT 5
TGGT 2
```
