# kcollections 3.x

Branch `v3.0-dev` — breaking release after v2.x.

## Done in 3.0.0 / 3.1.0

- [x] Remove deprecated Python API (`write`/`read` as primary API, trie debug on containers)
- [x] Serialization **v2**: little-endian `kcollections-v2` (file version 2)
- [x] Slimmer `Kdict` bindings (scalar + `vector_*` only)
- [x] Python **3.10+** only
- [x] Drop vendored `uint256_t` — `PrefMask` in `inc/pref_mask.h`
- [x] Dropped `uint256_t` (`PrefMask`)
- [x] Factored C++ layout: per-kind symbols (`KcSetContainer`, `KcDictContainer<T>`, …) in one `_kcollections` module
- [x] `std::thread` / `std::mutex` parallel add (replaces `pthread` + `sem_open`)
- [x] `add_seq` accepts `str`, `bytes`, `bytearray`, and buffer objects
- [x] Text `export_kmers` / `import_kmers` helpers
- [x] CI without Boost; benchmark smoke job
- [x] `CONTRIBUTING.md`

## Planned

- [ ] Fix `parallel_add` after `pthread` → `std::thread` migration (test skipped)
- [ ] Read v1 archives for one-release migration helper
- [ ] PyPI wheels + Bioconda
- [ ] Remove vendored `libs/pybind11-2.4.3/` from git history

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
