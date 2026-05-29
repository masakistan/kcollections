# Conda / Bioconda packaging

This recipe is a template for submitting **kcollections** to [Bioconda](https://bioconda.github.io/).

## Before submitting

1. Tag a release on GitHub, e.g. `v3.3.0`.
2. Update `sha256` in `meta.yaml`:

   ```bash
   url="https://github.com/masakistan/kcollections/archive/refs/tags/v3.3.0.tar.gz"
   curl -sL "$url" | sha256sum
   ```

3. Copy `meta.yaml` into a new folder under [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes) and open a PR there (Bioconda does not build from this repo directly).

## Local smoke test (optional)

```bash
conda install -c conda-forge conda-build
conda build conda-recipe
```
