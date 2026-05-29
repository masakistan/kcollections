# Usage guide

## When to use kcollections

| Goal | Use kcollections? | Alternative |
|------|-------------------|-------------|
| In-memory Python `set`/`dict` for k-mers, RAM-limited | **Yes** | Built-in `set` if &lt; ~1M k-mers |
| Attach metadata per k-mer (`Kdict`) | **Yes** | Rare elsewhere in one API |
| Prototype a method in Python | **Yes** | — |
| Billion-scale on-disk index | No | [KMC](https://github.com/KmerTools/KMC), Jellyfish |
| Sketching / MinHash | No | sourmash, mash |
| Production read mapper | No | BWA, minimap2, etc. |

## Recommended API (2.1+)

```python
from kcollections import Kset, Kdict, Kcounter

# Set of k-mers
ks = Kset(27)
ks.add_seq(dna)
assert "ACGT..." in ks

# Save / load (preferred names)
ks.save("index.kc")
ks2 = Kset.from_file("index.kc")

# Dict with values
kd = Kdict(int, 27)
kd["ACGT..."] = 42

# Counter
kc = Kcounter(27)
kc.add_seq(dna)
kc.most_common(10)
```

## Parallel insertion

Thread count must be a **power of 2** (implementation detail).

```python
ks = Kset(27)
ks.parallel_add_init(16)
ks.parallel_add_seq(dna)
ks.parallel_add_join()
```

## Persistence notes (production)

- Use `save()` / `load()` or `from_file()` — aliases for `write()` / `read()`.
- Indexes use the native **`kcollections-v1`** format (no Boost). Pin versions; 2.0/2.1 Boost archives are incompatible with 2.2+.
- Pin `kcollections==X.Y.Z` in production requirements.
- For long-term archival across machines, export k-mers to a simple text or FASTA-derived format and rebuild the index.

## Research workflows

1. **Index a reference region** — `Kdict(list, k)` + `add_seq` with locus lists (see `examples/dict_loci.py`).
2. **Compare abundance** — `Kcounter` + `most_common()`.
3. **Membership at scale in RAM** — `Kset` + `parallel_add_seq` on whole chromosomes.
4. **Debug trie structure** — `kcollections.debug.inspect(ks)` (internal; may change).

## Loading a saved `Kdict`

`k` and value types are restored from the file; pass a placeholder type when constructing:

```python
kd = Kdict(int, 0)
kd.load("index.kc")
```

Or:

```python
from kcollections import kdict_from_file
kd = kdict_from_file(int, "index.kc")
```
