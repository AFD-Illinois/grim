#include "grim.h"
#include <yaml-cpp/yaml.h>
#include <arrayfire.h>


int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


  DM dm;

  DMBoundaryType boundaryX1 = DM_BOUNDARY_GHOSTED;
  DMBoundaryType boundaryX2 = DM_BOUNDARY_GHOSTED;

  int numX1 = 4;
  int numX2 = 4;
  int numVar = 2;

  DMDACreate2d(PETSC_COMM_WORLD, 
               boundaryX1, boundaryX2,
               DMDA_STENCIL_BOX,
               numX1, numX2,
               PETSC_DECIDE, PETSC_DECIDE,
               numVar, 0,
               PETSC_NULL, PETSC_NULL,
               &dm
              );

  Vec vec;
  DMCreateGlobalVector(dm, &vec);

  double ***host_array;
  DMDAVecGetArrayDOF(dm, vec, &host_array);

  for (int j=0; j<numX2; j++)
  {
    for (int i=0; i<numX1; i++)
    {
      for (int var=0; var<numVar; var++)
      {
        host_array[j][i][var] = var;
      }
    }
  }

  af::info();
  af::array a(numVar, numX1, numX2, &(host_array[0][0][0]));
  af_print(af::shift(a(af::span, 0, 0), -1) );

  DMDAVecRestoreArrayDOF(dm, vec, &host_array);

  PetscFinalize();  
  return(0);
}
