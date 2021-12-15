#!/bin/bash
#BSUB -P AST171
#BSUB -W 00:10
#BSUB -J GRIM
# debug or batch
#BSUB -q debug

#BSUB -nnodes 1
NNODES=1

# Make bsub behave like sbatch, moving to CWD
# But don't cd anywhere under $HOME, that would be bad
if [[ "$(realpath $LS_SUBCWD)" !=  *"home"* ]]; then
  cd $LS_SUBCWD
fi

# Stuff for posterity
date
echo "Job run on nodes:"
jsrun -n $NNODES -r 1 hostname

module load gcc cuda/11.4 petsc boost nsight-systems

# Arrayfire configuration:
# Limit kernel lengths for performance
export AF_CUDA_MAX_JIT_LEN=30
# Cache to a directory the compute nodes can actually *see*
export AF_JIT_KERNEL_CACHE_DIRECTORY=$(realpath $PWD)/.arrayfire
# Debugging output
#export AF_TRACE=all

# Tell OpenMP explicitly what it's dealing with
# This also ensures that if we set e.g. NTHREADS=6,
# they will be spread 1/core
# Remember to change this if using 7 cores!
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=24

# In order to run anything on the compute nodes, we must specify
# how we wish to slice them up.  This is done with arguments to
# 'jsrun', which takes the place of 'mpirun' on Summit
# -n # resource sets (== MPI tasks for us)
# -r # rs per node
# -a # MPI task per rs
# -g # GPU per rs
# -c # CPU cores (physical) per rs (total of 42=(22-1)*2, but 7 is such an awkward number)
# -b binding strat *within* a resource set
# -l latency optimization in choosing cores (CPU-CPU, GPU-CPU, CPU-MEM)

# The "smpiargs" argument is used to tell Spectrum MPI we're GPU-aware
jsrun --smpiargs="-gpu" -n $(( $NNODES * 6 )) -r 6 -a 1 -g 1 -d packed -c 6 -b packed:6 -l GPU-CPU \
      nsys profile --capture-range=cudaProfilerApi -o grim_%q{OMPI_COMM_WORLD_RANK} ~/grim/grim