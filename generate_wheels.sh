#!/bin/bash
set -e -x

sudo yum install centos-release-scl;
sudo yum install llvm-toolset-7;
scl enable llvm-toolset-7 bash;
export CC=clang-7;
export CXX=clang++-7;

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install .
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install kcollections --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/nosetests" kcollections)
done
