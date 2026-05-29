# Contributing to kcollections

## Development setup

Requirements: Python 3.10+, CMake 3.18+, a C++17 compiler.

```bash
git clone https://github.com/masakistan/kcollections.git
cd kcollections
pip install -e ".[dev]"
pytest tests -q
```

No Boost or vendored pybind11 is required; the build uses `pybind11` from PyPI via scikit-build-core.

### C++ layout

The trie core is compiled once per kind with distinct type names (see `inc/kc/kset_core.h`, `kdict_core.h`, `kcounter_core.h`). Python bindings live in `kcollections/src/bindings_*.cc` and are linked into a single `_kcollections` extension.

## Build options

| CMake option | Default | Purpose |
|--------------|---------|---------|
| `KCOLLECTIONS_OPENMP` | OFF | Optional OpenMP (experimental) |
| `KCOLLECTIONS_BUILD_CPP_TESTS` | OFF | Legacy C++ test binaries |

## Pull requests

- Target `master` for 3.x maintenance and releases.
- Run `pytest tests -q` before opening a PR.
- Keep serialization format changes documented in `CHANGELOG.md` and `MIGRATION.md`.
- Prefer focused diffs; avoid drive-by refactors.

## Code style

- Match existing naming in `inc/` and `kcollections/`.
- C++: `-Wall -Wextra`, C++17.
- Python: type hints on public APIs in `kcollections/__init__.py`.

## Reporting issues

Include Python version, OS, k-mer length `k`, and a minimal reproducer (ideally without BioPython unless the bug is in optional `bio` extras).
