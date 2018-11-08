#!/usr/bin/env bash

set -e -x

# Install a system package required by our library
yum install zlib-devel -y
curl -o cmake.tar.gz https://cmake.org/files/v3.12/cmake-3.12.4.tar.gz
tar xf cmake.tar.gz
cd cmake-3.12.4/
./bootstrap --system-curl
gmake install

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/python" setup.py bdist_wheel sdist
done
cp /io/dist/*.tar.gz /io/wheelhouse

# Bundle external shared libraries into the wheels
for whl in /io/dist/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install kcollections --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/python" -c 'from kcollections import Kset, Kdict')
done