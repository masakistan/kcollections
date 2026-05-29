# kcollections 3.x

Maintenance on `master`.

## Done (3.0 – 3.3)

- [x] Serialization v2; unified `_kcollections` module; `PrefMask`; Python 3.10+
- [x] `parallel_add` + work-queue drain fix; `debug.inspect()`; text export/import
- [x] Read **v1** archives (LE) and `python -m kcollections migrate` helper (3.3)
- [x] PyPI publish workflow on version tags (3.3)
- [x] Remove vendored `libs/pybind11-2.4.3/` (3.3.1)
- [x] Bioconda recipe template under `conda-recipe/` (3.3.1)
- [x] Native C++ smoke test (`KCOLLECTIONS_BUILD_CPP_TESTS`) + CI job (3.3.1)

## Planned

- [ ] Submit recipe to [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) (update `sha256` per release)
- [ ] Optional: rewrite git history to drop old `libs/pybind11` blobs (large clone size)

## Migration

**2.2 v1 → 3.x v2** (little-endian):

```bash
python -m kcollections migrate old.kc new.kc
```

Boost-era 2.0/2.1 files are not supported — rebuild indexes.
