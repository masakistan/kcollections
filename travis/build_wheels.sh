#!/usr/bin/env bash

set -e -x

# Install a system package required by our library
yum install zlib-devel -y
yum install python-devel -y

echo "installing boost"
curl -L -o boost.tar.gz https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.gz/download
tar xf boost.tar.gz
cd boost_1_72_0

./bootstrap.sh --with-libraries=serialization --prefix=/usr && \
    ./b2 stage -j1 threading=multi link=shared

./b2 install threading=multi link=shared && \
    ln -svf detail/sha1.hpp /usr/include/boost/uuid/sha1.hpp
echo "done installing boost"

curl -o cmake.tar.gz https://cmake.org/files/v3.12/cmake-3.12.4.tar.gz
tar xf cmake.tar.gz
cd cmake-3.12.4/
pwd
./bootstrap --system-curl
gmake install

cd /io
echo "checking /io"
pwd
ls

export PYBIN=/opt/python/cp"$1"-cp"$1"mu/bin

# Compile wheels
"${PYBIN}/python" /io/setup.py bdist_wheel sdist
cp /io/dist/*.tar.gz /io/wheelhouse
echo "checking /io/wheelhouse"
ls /io/wheelhouse

# Bundle external shared libraries into the wheels
for whl in /io/dist/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

# Install packages and test
"${PYBIN}/pip" install kcollections --no-index -f /io/wheelhouse
(cd "$HOME"; "${PYBIN}/python" -c 'from kcollections import Kset, Kcounter, Kdict_int, Kdict_float, Kdict_string')
