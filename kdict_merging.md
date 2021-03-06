# Kdict

## Declaring a Kdict

Because the C++ bindings are compiled in `kcollections` the value type must be specified when using `Kdict`.
See the last section of this page for implemented types.
In the following example we create a `Kdict` that will store 4-mers as keys and ints as values:

```python
import kcollections
kd = kcollections.Kdict(int, 4)
```

The following example initializes a `Kdict` with values that are lists of ints.

```python
import kcollections
kd = kcollections.Kdict((list, int), 4)
```


## The Merging Function
`Kdict` provides a merginging function as a means to reconcile two values for a kmer that are being added to the data structure.
The default merging function is:

```python
def merge_func(prev_val, new_val):
    return new_val
```

This simply replaces the old value for a kmer with the new value.

This is not interesting but given a different merging function more interesting information can be stored in a `Kdict` as seen in the following examples.


## Examples

### Default Merging Function
This example uses the default mergining function which simply overwrites the previous value.
It stores the last index of each kmer in the dna string.

```python
import kcollections
dna = 'ACTGGTACTG'

kd = kcollections.Kdict(int, 4)

kd.parallel_add_init(4)
kd.parallel_add_seq(dna, [i for i in range(len(dna))])
kd.parallel_add_join()

for kmer, val in kd.iteritems():
    print(kmer, val)
```

Output:

```bash
GGTA 3
GTAC 4
CTGG 1
ACTG 6
TACT 5
TGGT 2
```

### Kmer Counting
Though we provide `Kcounter`, kmer counting could be accomplished using `Kdict` with no difference in performance.

``` python
import kcollections
dna = 'ACTGGTACTG'

kd = kcollections.Kdict(int, 4)
kd.set_merge_func(lambda prev_val, new_val: prev_val + new_val)

kd.parallel_add_init(4)
kd.parallel_add_seq(dna, [1 for _ in range(len(dna))])
kd.parallel_add_join()

for kmer, val in kd.iteritems():
    print(kmer, val)
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
The following example stores all indices of each kmer occurrence in a DNA string.

``` python
import kcollections
dna = 'ACTGGTACTG'

def merge_func(prev_val, new_val):
    prev_val.append(new_val[0])
    return prev_val

kd = kcollections.Kdict((list, int), 4)
kd.set_merge_func(merge_func)
kd.parallel_add_init(4)
kd.parallel_add_seq(dna,[[i] for i in range(len(dna))])
kd.parallel_add_join()

for kmer, val in kd.items():
    print(kmer, val)
```

Output:

``` bash
GGTA [3]
GTAC [4]
CTGG [1]
ACTG [0, 6]
TACT [5]
TGGT [2]
```

## PyObject Values and the GIL

`pybind11` provides automatic type conversion for many types and we have implemented several cases but there may be instances where you want to use a custom object.
While it is possible to create a `Kdict` with an `py::object` value, we suggest that you implement a custom class in C++.
Using `object` as the value type will cause `Kdict` to acquire the Global Interpreter Lock (GIL) when performing operations which negates all speedup from parallelization.
The use of the GIL makes parallel operatoins (essentially) sequential removing any benefit of parallel operation.

## List of Available Value Types
A type of `(list, int)` will create a list of ints as the value

- `int`
- `float`
- `str`
- `bool`
- `(list, int)`
- `(list, float)`
- `(list, str)`
- `(list, bool)`
