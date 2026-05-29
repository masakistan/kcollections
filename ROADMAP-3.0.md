# kcollections 3.0

Branch `v3.0-dev` — breaking release after v2.x.

## Done in 3.0.0 (this branch)

- [x] Remove deprecated Python API (`write`/`read`, `iteritems`, `Kdict_*`, trie debug on containers)
- [x] Serialization **v2**: little-endian `kcollections-v2` (file version 2)
- [x] Slimmer `Kdict` bindings (scalar + `vector_*` only; no `set_`/`list_`/`list_list`)
- [x] Python **3.10+** only

## Planned (3.0.x / 3.1)

- [ ] Single extension module `_kcollections`
- [ ] Read v1 archives for one-release migration helper
- [ ] PyPI wheels + Bioconda
- [ ] Delete vendored `libs/pybind11-2.4.3/`
- [ ] Optional `numpy` buffer `add_seq`

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
