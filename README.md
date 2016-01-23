grim
====

General Relativistic Implicit Magnetohydrodynamics

Instructions to compile PETSc on UIUC campus cluster:

1) Get L1 cacheline size:

$ getconf LEVEL1_DCACHE_LINESIZE

2) Get L1 cache size

either

$ grep . /sys/devices/system/cpu/cpu0/cache/index*/*

or 

$ lstopo-no-graphics

3) To get L1 cache associativity:

$ grep . /sys/devices/system/cpu/cpu0/cache/index*/*

4) Compiling petsc
./configure --prefix=/home/manic/petsc_optimized/ --with-debugging=0 COPTFLAGS='-O3 -qopt-report=5  -qopt-report-phase=vec -xhost' CXXOPTFLAGS='-O3 -qopt-report=5  -qopt-report-phase=vec -xhost' --with-hdf5=1 --with-clean=1 --with-mpi-dir=/usr/local/mpi/openmpi-1.8.4-intel-15.0/ --with-memalign=64 --known-level1-dcache-size=32768 --known-level1-dcache-linesize=64 --known-level1-dcache-assoc=8
