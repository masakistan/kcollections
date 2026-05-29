# Legacy material

This directory is **not maintained** and exists only for historical reference.

- `scripts/` — Python 2-era benchmarks and experiments from the original development period. Use [`examples/`](../examples/) and [`tests/`](../tests/) instead.

The build uses **pybind11 from PyPI** (via scikit-build-core). Vendored `libs/pybind11-2.4.3/` was removed in 3.3.1; `libs/uint256_t/` remains only as historical reference (replaced by `inc/pref_mask.h`).
