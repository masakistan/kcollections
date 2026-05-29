# kcollections 3.x

Maintenance on `master` after the 3.2 release.

## Done (3.0 – 3.2)

- [x] Remove deprecated Python API (`write`/`read` as primary API, trie debug on containers)
- [x] Serialization **v2**: little-endian `kcollections-v2` (file version 2)
- [x] Slimmer `Kdict` bindings (scalar + `vector_*` only)
- [x] Python **3.10+** only
- [x] `PrefMask` replaces vendored `uint256_t`
- [x] Factored C++ layout + single `_kcollections` module (3.2)
- [x] `std::thread` / `std::mutex` parallel add; `parallel_add` test fixed (3.2.1)
- [x] `add_seq` accepts `str`, `bytes`, `bytearray`, and buffer objects
- [x] Text `export_kmers` / `import_kmers` helpers
- [x] CI without Boost; benchmark smoke job
- [x] `CONTRIBUTING.md`
- [x] `kcollections.debug.inspect()` via internal `_get_*` bindings (3.2.1)

## Planned

- [ ] Read v1 archives for one-release migration helper
- [ ] PyPI wheels + Bioconda (cibuildwheel on tag push)
- [ ] Remove vendored `libs/pybind11-2.4.3/` from git history
- [ ] Re-enable C++ unit tests for factored layout

## Migration 2.2 → 3.0

1. Upgrade Python to 3.10+.
2. Replace deprecated calls (already warned in 2.1–2.2).
3. **Re-save all indexes** — v1 (`kcollections-v1`) files are not readable in 3.0.

```python
from kcollections import Kset

ks = Kset(27)
# ... populate ...
ks.save("index.kc")  # writes v2 format
```
