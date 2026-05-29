# Release checklist

## 1. Version bump

Update `pyproject.toml` and `kcollections/__init__.py`, add a section to `CHANGELOG.md`.

## 2. Tag and push

```bash
git tag -a v3.3.1 -m "kcollections 3.3.1"
git push origin v3.3.1
```

Pushing a `v*` tag triggers:

- **Publish** workflow — cibuildwheel + PyPI (requires [trusted publisher](https://docs.pypi.org/trusted-publishers/) on PyPI for this repo)
- **Release** workflow — GitHub release with generated notes

## 3. Bioconda

After the tag exists:

```bash
./scripts/conda_sha256.sh 3.3.1
```

Copy `conda-recipe/meta.yaml` into [bioconda-recipes/recipes/kcollections/](https://github.com/bioconda/bioconda-recipes/tree/master/recipes) with the updated `sha256` and open a PR.

## 4. Verify

```bash
pip install kcollections==3.3.1
pytest tests -q
python -m kcollections probe tests/fixtures/v1_kset_seq.kc
```
