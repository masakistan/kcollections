# Migration guide

## Goals in 2.1+

- **Minimal maintenance** — small public API, examples + pytest, legacy code quarantined under `legacy/`
- **Ease of use** — `save`/`load`/`from_file`, one import path
- **Production** — pin versions; understand serialization limits (below)
- **Research** — `Kdict` with typed values, `examples/`, optional `kcollections.debug`

## From 1.x / Python 2

| Old | New |
|-----|-----|
| Python 2.7 | Python 3.8+ |
| `pip install kcollections` (old wheels) | Build from source or wait for 2.1+ PyPI wheels |
| `print kmer` | `print(kmer)` |
| `kd.iteritems()` | `kd.items()` |
| `Kdict_int(27)` | `Kdict(int, 27)` |
| `ks.write(path)` | `ks.save(path)` (preferred) |
| `scripts/*.py` | `examples/*.py` |

## Deprecated in 2.1 (removed in 3.0)

| API | Use instead |
|-----|-------------|
| `write()` / `read()` | `save()` / `load()` or `from_file()` |
| `iteritems()` | `items()` |
| `Kdict_int`, `Kdict_float`, … | `Kdict(int, k)`, etc. |
| `get_root()`, `get_uc_kmer()`, … on containers | `kcollections.debug.inspect(obj)` |
| `Kdict((list, list), k)` nested containers | Scalar or `Kdict(list, k)` only |

Enable warnings in notebooks:

```python
import warnings
warnings.simplefilter("always", DeprecationWarning)
```

## When not to use kcollections

Use **KMC**, **Jellyfish**, or **Cuttlefish** for huge on-disk k-mer databases in pipelines. Use **kcollections** for in-memory Python prototypes and research tools where `dict`/`set` semantics matter.

## Serialization compatibility

**2.2+** uses native format `kcollections-v1` (magic `KCOL`, no Boost). Files written by **2.0/2.1** (Boost) must be re-exported or rebuilt.

Safe when:

- Same `kcollections` major version (2.2.x)
- Same container type (`Kset` vs `Kdict` vs `Kcounter`)
- Same `Kdict` value type (e.g. `int` vs `str`)

Endianness follows the platform that wrote the file. For archival, export k-mers to text and rebuild.
