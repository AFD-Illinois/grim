#!/bin/bash

# Compile grim, assuming we're on Summit and the dependencies are available
# USE:
# ./make.sh [clean]

NPROC=
if [[ $(hostname) == "login"* ]]; then
  NPROC=16
fi

if [ ! -d build ]; then
  echo "Run \"./make_summit.sh clean\" first!"
  exit
fi

module load gcc cuda/11.4 petsc boost python

module list

if [[ "$*" == *"clean"* ]]; then
  rm -rf build
  mkdir build
  cd build
  cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
        -DPETSC_EXECUTABLE_RUNS=1 \
        -DArrayFire_ROOT_DIR=/gpfs/alpine/proj-shared/ast171/libs/arrayfire-cuda114 \
        ../src
  cd ..
fi

cd build
make -j$NPROC
cp grim ..
