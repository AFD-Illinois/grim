#include "grim.hpp"
#include "params.cpp"

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscPrintf(PETSC_COMM_WORLD, 
              "#### System info ####\n"
             );
  if (rank==0)
  {
    af::info();
  };

  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    timeStepper ts;

    params::Time = 0.;
    while(params::Time < params::finalTime)
    {
    	ts.timeStep();
    	params::Time+=params::dt;
     }
  }

  PetscFinalize();  
  return(0);
}
