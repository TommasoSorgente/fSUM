#!/bin/bash

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BUILDDIR=build
SYSBUILDDIR=$BUILDDIR-$OSTYPE

echo $SCRIPTDIR

cd $SCRIPTDIR/segmentation
mkdir -p $SYSBUILDDIR
cd $SYSBUILDDIR
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --parallel 8
cd ..
mv $SYSBUILDDIR $BUILDDIR

cd $SCRIPTDIR/merge_subdomains
mkdir -p $SYSBUILDDIR
cd $SYSBUILDDIR
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel 8
cd ..
mv $SYSBUILDDIR $BUILDDIR

cd $SCRIPTDIR

