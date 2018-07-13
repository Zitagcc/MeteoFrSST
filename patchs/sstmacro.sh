#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Usage: $0 release|debug install_dir"
  exit
fi

TOP=`pwd`
mkdir $TOP/b
cd $TOP/b

if [ "$1" == "release" ]; then
  ../configure --prefix=$2 CC=mpicc CXX=mpicxx
else
  ../configure --prefix=$2 CC=mpicc CXX=mpicxx CFLAGS="-g -O0" CXXFLAGS="-g -O0"
fi

make -j 16 install

cd $TOP
