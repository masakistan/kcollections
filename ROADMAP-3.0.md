# kcollections 3.x

Maintenance on `master`.

## Done (3.0 – 3.3)

- [x] Serialization v2; unified `_kcollections` module; `PrefMask`; Python 3.10+
- [x] `parallel_add` + work-queue drain fix; `debug.inspect()`; text export/import
- [x] Read **v1** archives (LE) and `python -m kcollections migrate` helper (3.3)
- [x] PyPI publish workflow on version tags (3.3)

## Planned

- [ ] Bioconda recipe
- [ ] Remove vendored `libs/pybind11-2.4.3/` from git history
- [ ] Re-enable C++ unit tests for factored layout

## Migration

**2.2 v1 → 3.x v2** (little-endian):

```bash
python -m kcollections migrate old.kc new.kc
# or in Python: load with Kset.from_file("old.kc"); save("new.kc")
```

Boost-era 2.0/2.1 files are not supported — rebuild indexes.
