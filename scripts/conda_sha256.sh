#!/usr/bin/env bash
# Print sha256 for a GitHub release tarball (Bioconda source.sha256).
set -euo pipefail

version="${1:?usage: $0 VERSION (e.g. 3.3.1)}"
url="https://github.com/masakistan/kcollections/archive/refs/tags/v${version}.tar.gz"
echo "url: ${url}"
curl -sL "${url}" | shasum -a 256
