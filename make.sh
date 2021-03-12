#!/bin/bash

# Download & compile dependencies of grim,
# then compile grim itself

# TODO
# Make YAML a submodule?
# Power9 will be a *nightmare*

ARRAYFIRE_VER=3.6.4
ARRAYFIRE_SCRIPT=ArrayFire-v${ARRAYFIRE_VER}_Linux_x86_64.sh
NPROC=12

if [[ "$*" == *"dep"* ]]; then
  rm ArrayFire*
  rm -rf external
  mkdir external

  # Arrayfire
  wget https://arrayfire.s3.amazonaws.com/$ARRAYFIRE_VER/$ARRAYFIRE_SCRIPT
  chmod +x $ARRAYFIRE_SCRIPT
  ./$ARRAYFIRE_SCRIPT --prefix=$PWD/external --include-subdir --skip-license

  # YAML
  cd external
  git clone https://github.com/jbeder/yaml-cpp
  cd yaml-cpp
  mkdir build
  cd build
  cmake -DYAML_BUILD_SHARED_LIBS=ON ..
  make -j$NPROC
  cd ../../..
fi

if [[ "$*" == *"clean"* ]]; then
  rm -rf build
  mkdir build
  cd build
  cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ../src
  #cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ../src
  cd ..
fi

cd build
make -j$NPROC
cp grim ..
