#!/bin/sh
echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

fn=metis-4.0.3.tar.gz
if [ ! -f $fn ]; then
  echo "### Downloading Metis:"
  if wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/${fn}
  then
    echo "### Metis: Download Successful.\n"
  else
    echo "### Metis: Download Failed.\n"
    exit 1
  fi
else
  echo "File found. No download."
fi

name=`basename ${fn} .tar.gz`
rm -rf ${name} src
tar -zxvf $fn
ln -s ./${name} ./src

#compile metis
cd src
sed -i  "s/\bCOPTIONS =/COPTIONS = -O3 -fPIC /g" Makefile.in
make -j$1


