#!/bin/bash

## This script should help to bootstrap PIPS.

# exit if a command fails
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

cmake -DDUMP=OFF -DELEMENTAL=OFF -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ALL=OFF -DBUILD_PIPS_NLP=ON -DBUILD_PIPS_DOC=ON -B. -H..
#cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ALL=OFF -DBUILD_PIPS_NLP=ON -DBUILD_PIPS_DOC=ON -B. -H..

# Build using 4 processes. 
make -j

