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

## 2.x → 3.0 (breaking)

| Removed in 3.0 | Use instead |
|----------------|-------------|
| `write()` / `read()` | `save()` / `load()` or `from_file()` |
| `iteritems()` | `items()` |
| `Kdict_int`, `Kdict_float`, … | `Kdict(int, k)`, etc. |
| `get_root()`, `get_uc_kmer()`, … on containers | `kcollections.debug.inspect(obj)` |
| `Kdict_set_*`, `Kdict_list_*`, `list_list` | `Kdict(int, k)` or `Kdict((list, int), k)` |
| Python 3.8–3.9 | Python 3.10+ |
| Archive v1 | Re-save indexes (v2 little-endian format) |

## When not to use kcollections

Use **KMC**, **Jellyfish**, or **Cuttlefish** for huge on-disk k-mer databases in pipelines. Use **kcollections** for in-memory Python prototypes and research tools where `dict`/`set` semantics matter.

## Serialization compatibility

**3.0+** uses `kcollections-v2` (magic `KCOL`, version 2, little-endian). Not compatible with 2.2 v1 archives or Boost-era 2.0/2.1 files — rebuild indexes.

Safe when:

- Same `kcollections` 3.x version
- Same container type (`Kset` vs `Kdict` vs `Kcounter`)
- Same `Kdict` value type (e.g. `int` vs `(list, int)`)
