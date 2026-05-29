# Changelog

## 3.0.0 (v3.0-dev)

### Added
- Serialization v2: little-endian `kcollections-v2` (portable across LE platforms)
- `ROADMAP-3.0.md`

### Changed
- Python 3.10+ required
- `Kdict` bindings: scalar + `(list, T)` / `vector_*` only

### Removed
- Deprecated `write`/`read`, `iteritems`, `Kdict_*` exports, trie methods on containers
- `Kdict` `set_`/`list_`/`list_list` C++ types

### Breaking
- v1 archive files from 2.2 must be re-saved
- See [MIGRATION.md](MIGRATION.md)

## 2.2.0

### Changed
- **Removed Boost dependency** — persistence uses native `kcollections-v1` binary format (`inc/kc_io.h`)
- Install requires only CMake + C++17 (no `libboost-serialization-dev`)

### Breaking
- On-disk files from 2.0/2.1 (Boost archives) are **not** readable in 2.2+; rebuild indexes or keep an older version to migrate

## 2.1.0

### Added
- `save()` / `load()` / `from_file()` as the preferred persistence API
- `kdict_from_file()` helper
- `kcollections.debug.inspect()` for trie debugging
- `examples/`, `docs/USAGE.md`, `MIGRATION.md`
- `py.typed` for type checkers
- Deprecation warnings with removal planned in 3.0

### Changed
- Public API narrowed: `Kset`, `Kdict`, `Kcounter` only in `__all__`
- Legacy `scripts/` moved to `legacy/scripts/` (unmaintained)
- README focused on research + light production use

### Deprecated
- `write()` / `read()`, `iteritems()`, `Kdict_*` exports, trie introspection on containers

## 2.0.0

### Added
- `pyproject.toml` with scikit-build-core for modern Python packaging
- GitHub Actions CI (Linux + macOS, Python 3.9–3.13)
- cibuildwheel configuration for manylinux and macOS wheels
- Pytest suite under `tests/`
- `Kcounter.most_common()` helper
- Default constructors on `Kset`/`Kdict`/`Kcounter` for `read()` workflows
- `thread_local` `CDEPTH` for safer concurrent serialization

### Changed
- Python 3.8+ only (dropped Python 2)
- C++17 standard; pybind11 2.11+ via build dependencies (no vendored pybind11 at build time)
- Portable `CHECK_KMER_LENGTH` (no GNU statement expressions)
- `Kdict` factory uses `getattr` instead of `eval()`
- README and examples updated for Python 3

### Fixed
- `Kset.intersection_update` / `symmetric_difference_update` now mutate in place
- `Kset.__or__` returns a new set instead of mutating in place
- `pop(key, default)` returns the default value, not a tuple
- CMake test targets link the correct libraries (`Kdict`, `Kcounter`)
- Removed exposed Travis CI token from README badge
- `applications/mapper.py` updated for current API

### Removed
- Travis CI configuration (replaced by GitHub Actions)
