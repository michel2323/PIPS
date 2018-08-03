#!/bin/bash

## This script should help to bootstrap PIPS.

# exit if a command fails
export CRAYPE_LINK_TYPE=dynamic
export LD_LIBRARY_PATH=/home/mschanen/git/PIPS/ThirdPartyLibs/Elemental/install/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/mschanen/git/PIPS/ThirdPartyLibs/Elemental/install/lib:$LD_LIBRARY_PATH

set -e
if [ -z $1 ]; then
  NUMTHREADS=4
else
  NUMTHREADS=$1
fi

rm -rf build
mkdir build
cd build
# CMAKE_BUILD_TYPE is RELEASE or DEBUG
# BUILD_ALL=OFF: Do not build all PIPS variants
# BUILD_PIPS_NLP: Build PIPS_NLP
# BUILD_PIPS_DOC: Build the documentation target (make doc)
# DUMP=ON: Dump 1st stage matrix

cmake -DELEMENTAL=OFF -DMATH_LIBS='-ldl' -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ALL=OFF -DBUILD_PIPS_NLP=ON -DBUILD_PIPS_DOC=ON -DDUMP=OFF -B. -H..

# Build using 4 processes. 
make -j 

