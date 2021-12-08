#!/bin/bash

# Compile grim, assuming we're on Summit and the dependencies are available
# USE:
# ./make.sh [clean]

module load gcc cuda petsc/3.15.0 boost

NPROC=
if [[ $(hostname) == "login"* ]]; then
  NPROC=16
fi

if [[ "$*" == *"clean"* ]]; then
  rm -rf build
  mkdir build
  cd build
  cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
        -DPETSC_EXECUTABLE_RUNS=1 \
        -DArrayFire_ROOT_DIR=/gpfs/alpine/proj-shared/ast171/libs/arrayfire \
        ../src
  cd ..
fi

cd build
make -j$NPROC
cp grim ..
