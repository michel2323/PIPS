#!/bin/sh
fn=Cbc-2.9.8.tgz

echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

if [ ! -f $fn ]; then
  echo "### Downloading Cbc:"
  if wget http://www.coin-or.org/download/source/Cbc/${fn}
  then
    echo "### ${fn}: Download Successful.\n"
  else
    echo "### ${fn}: Download Failed.\n"
    exit 1
  fi
else
  echo "File found. No download."
fi
rm -rf Cbc-2.9.8 src
name=`basename ${fn} .tgz`
tar -zxf $fn
git apply ${name}.patch
ln -s ./${name} ./src

cd src
#./configure --enable-static --prefix=`pwd`
#make -j4 install
./configure --prefix=`pwd` --disable-shared --enable-static=yes --host x86_64 CFLAGS="-fPIC -O3" CXXFLAGS="-fPIC -O3"
make -j install

