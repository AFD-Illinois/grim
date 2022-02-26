#include "grim.hpp"
#include "params.hpp"

//#include "cuda_profiler_api.h"

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);

  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  af::setDevice(world_rank%params::numDevices);

  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    timeStepper ts(params::N1, params::N2, params::N3,
                   params::dim, vars::dof, params::numGhost,
                   params::Time, params::InitialDt,
                   params::boundaryLeft,  params::boundaryRight,
                   params::boundaryTop,   params::boundaryBottom,
                   params::boundaryFront, params::boundaryBack,
                   params::metric, params::blackHoleSpin, params::hSlope,
                   params::X1Start, params::X1End,
                   params::X2Start, params::X2End,
                   params::X3Start, params::X3End
                  );

    PetscPrintf(PETSC_COMM_WORLD, "  Generating compute kernels...\n\n");
    int numReads, numWrites;
    ts.timeStep(numReads, numWrites);

    af::sync();

    PetscPrintf(PETSC_COMM_WORLD, "\n  Kernel compilation complete\n");

    //cudaProfilerStart();

    int n=0;
    int StopRunning = 0;
    while(ts.time<params::finalTime && StopRunning==0)
    {
      n++;
      PetscPrintf(PETSC_COMM_WORLD, "\n|----Time step %d----|\n", n);
      ts.timeStep(numReads, numWrites);
      //Checkpoint if running out of time
      StopRunning = ts.CheckWallClockTermination();
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n===Program execution complete===\n\n");
    if(StopRunning)
      PetscPrintf(PETSC_COMM_WORLD, "\n Termination reason: WallClock\n");
    else
      PetscPrintf(PETSC_COMM_WORLD, "\n Termination reason: Final Time\n");
  }
  PetscFinalize();
  return(0);
}
