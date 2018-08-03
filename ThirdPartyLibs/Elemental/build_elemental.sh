#!/bin/sh

set -e
export CRAYPE_LINK_TYPE=dynamic

rm -rf Elemental-0.87.7 install
wget https://github.com/elemental/Elemental/archive/v0.87.7.tar.gz -O v0.87.7.tar.gz
tar xzf v0.87.7.tar.gz
cd Elemental-0.87.7
mkdir build
cd build
cmake .. -DMATH_LIBS="-ldl" -DCMAKE_INSTALL_PREFIX=$PWD/../../install
#cmake .. -DMATH_LIBS="-Wl,--no-as-needed -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -lgfortran" -DCMAKE_INSTALL_PREFIX=$PWD/../../install
make -j
make install

