#!/bin/sh

set -e

wget https://github.com/elemental/Elemental/archive/v0.87.7.tar.gz -O v0.87.7.tar.gz
tar xzf v0.87.7.tar.gz
cd Elemental-0.87.7
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/../../
make -j$1
make install

