#!/bin/bash

# Download & compile dependencies of grim,
# then compile grim itself
# USE:
# ./make.sh [dep] [clean]

# Needs local MPI pointed to with MPI_DIR.  Most modules do this.

# TODO
# Tune PETSc compile
# Make YAML a submodule?
# Power9 will be a *nightmare*

ARRAYFIRE_VER=3.6.4
ARRAYFIRE_SCRIPT=ArrayFire-v${ARRAYFIRE_VER}_Linux_x86_64.sh
PETSC_VER=3.15.4
NPROC=56
USE_INTEL=0 # Bash comes with some tradeoffs.  Booleans are one of them.

if [[ $(hostname -f) == "bh"* ]]; then
  module load gnu mpich phdf5
fi

if [[ "$*" == *"arrayfire"* ]]; then
  mkdir -p external
  cd external

  rm -rf arrayfire-master
  git clone github:arrayfire/arrayfire arrayfire-master
  cd arrayfire-master

  mkdir build
  cd build

  cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
        -DAF_BUILD_CPU=OFF -DAF_BUILD_OPENCL=OFF -DAF_BUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/Code/grim/external/arrayfire ..
        #-DBOOST_ROOT=$HOME/grim/external/boost_1_76_0 -DBOOST_INCLUDEDIR=$HOME/grim/external/boost_1_76_0/boost \
        #-DBoost_NO_SYSTEM_PATHS=ON ..
  
  make -j$NPROC
  
  make install
  cd ../..
fi

if [[ "$*" == *"petsc"* ]]; then
  cd external
  rm -rf petsc-${PETSC_VER} petsc-${PETSC_VER}.tar.gz
  wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VER}.tar.gz
  tar xf petsc-${PETSC_VER}.tar.gz
  cd petsc-${PETSC_VER}

  if [[ $USE_INTEL != "0" ]]; then
    internal_COPTFLAGS='-O3 -qopt-report=5 -qopt-report-phase=vec -xhost'
    internal_CXXOPTFLAGS='-O3 -qopt-report=5 -qopt-report-phase=vec -xhost'
  else
    internal_COPTFLAGS='-O3 -march=native -mtune=native -funroll-loops'
    internal_CXXOPTFLAGS='-O3 -march=native -mtune=native -funroll-loops'
  fi

  ./configure --prefix=$PWD/../petsc --with-debugging=0 \
  COPTFLAGS=${internal_COPTFLAGS} CXXOPTFLAGS=${internal_CXXOPTFLAGS} \
  --with-hdf5=1 --with-clean=1 --with-mpi-dir=$(dirname $(which mpicc))/.. \
  --download-fblaslapack
  #--with-memalign=64 --known-level1-dcache-size=32768 --known-level1-dcache-linesize=64 \
  #--known-level1-dcache-assoc=8
  #python config/examples/arch-olcf-summit-opt.py

  make -j$NPROC
  make install

  cd ../..
fi



if [[ "$*" == *"clean"* ]]; then
  rm -rf build
  mkdir build
  cd build
  if [[ $USE_INTEL != "0" ]]; then
    cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ../src
  else
    cmake -DPETSC_EXECUTABLE_RUNS=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ../src
  fi
  cd ..
fi

cd build
make -j$NPROC
cp grim ..

echo "On BH, run the following to load MKL correctly"
echo ". /opt/intel/oneapi/setvars.sh"
