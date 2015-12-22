#include "grid.hpp"

grid::grid(int numVars, int numGhost)
{
  this->numVars = numVars;
  this->numGhost = numGhost;

  /* Implementations for MIRROR, OUTFLOW in boundary.cpp and DIRICHLET in
   * problem.cpp */
  boundaryLeft  = DM_BOUNDARY_GHOSTED; boundaryRight  = DM_BOUNDARY_GHOSTED;
  boundaryTop   = DM_BOUNDARY_GHOSTED; boundaryBottom = DM_BOUNDARY_GHOSTED;
  boundaryFront = DM_BOUNDARY_GHOSTED; boundaryBack   = DM_BOUNDARY_GHOSTED;

  if (   params::boundaryLeft  == boundaries::PERIODIC 
      || params::boundaryRight == boundaries::PERIODIC
     )
  {
    boundaryLeft  = DM_BOUNDARY_PERIODIC;
    boundaryRight = DM_BOUNDARY_PERIODIC;
  }

  if (   params::boundaryTop    == boundaries::PERIODIC 
      || params::boundaryBottom == boundaries::PERIODIC
     )
  {
    boundaryTop    = DM_BOUNDARY_PERIODIC;
    boundaryBottom = DM_BOUNDARY_PERIODIC;
  }

  if (   params::boundaryFront == boundaries::PERIODIC 
      || params::boundaryBack  == boundaries::PERIODIC
     )
  {
    boundaryFront = DM_BOUNDARY_PERIODIC;
    boundaryBack  = DM_BOUNDARY_PERIODIC;
  }

  switch (params::dim)
  {
    case 1:
      DMDACreate1d(PETSC_COMM_WORLD, boundaryLeft, 
                   params::N1, numVars, 0, NULL,
                   &dm
                  );

      DMDACreate1d(PETSC_COMM_WORLD, boundaryLeft, 
                   params::N1, numVars, numGhost, NULL,
                   &dmGhost
                  );

      break;
  
    case 2:
      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, 0,
                   PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost,
                   PETSC_NULL, PETSC_NULL,
                   &dmGhost
                  );
      break;

    case 3:
      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2, params::N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, 0, 
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2, params::N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost, 
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dmGhost
                  );
      break;
  }

  DMCreateGlobalVector(dm, &globalVec);
  DMCreateLocalVector(dmGhost, &localVec);
  VecSet(globalVec, 0);
  VecSet(localVec, 0);

  if (!gridParams::haveGridParamsBeenSet)
  {
    DMDAGetCorners
      (dm,
       &gridParams::iLocalStart,
       &gridParams::jLocalStart,
       &gridParams::kLocalStart,
       &gridParams::N1Local,
       &gridParams::N2Local,
       &gridParams::N3Local
      );

    gridParams::iLocalEnd = gridParams::iLocalStart + gridParams::N1Local;
    gridParams::jLocalEnd = gridParams::jLocalStart + gridParams::N2Local;
    gridParams::kLocalEnd = gridParams::kLocalStart + gridParams::N3Local;

    gridParams::dX1 = (params::X1End - params::X1Start)/params::N1;
    gridParams::dX2 = (params::X2End - params::X2Start)/params::N2;
    gridParams::dX3 = (params::X3End - params::X3Start)/params::N3;
  }

  vars = new array[numVars];
  array varsCopiedFromVec;

  if (numGhost > 0)
  {
    DMGlobalToLocalBegin(dmGhost, globalVec, INSERT_VALUES, localVec);
    DMGlobalToLocalEnd(dmGhost, globalVec, INSERT_VALUES, localVec);

    double *pointerToLocalVec;
    VecGetArray(localVec, &pointerToLocalVec);

    switch (params::dim)
    {
      case 1:
        varsCopiedFromVec = 
          array(numVars,
                gridParams::N1Local + 2*numGhost, 
                gridParams::N2Local,
                gridParams::N3Local,
                pointerToLocalVec
               );

        break;

      case 2:
        varsCopiedFromVec = 
          array(numVars,
                gridParams::N1Local + 2*numGhost, 
                gridParams::N2Local + 2*numGhost,
                gridParams::N3Local,
                pointerToLocalVec
               );

        break;

      case 3:
        varsCopiedFromVec = 
          array(numVars,
                gridParams::N1Local + 2*numGhost, 
                gridParams::N2Local + 2*numGhost,
                gridParams::N3Local + 2*numGhost,
                pointerToLocalVec
               );

        break;
      }

    VecRestoreArray(localVec, &pointerToLocalVec);
  }
  else
  {
    double *pointerToGlobalVec;
    VecGetArray(globalVec, &pointerToGlobalVec);

    varsCopiedFromVec = 
      array(numVars,
            gridParams::N1Local, gridParams::N2Local, gridParams::N3Local,
            pointerToGlobalVec
           );

    VecRestoreArray(globalVec, &pointerToGlobalVec);
  }

  array varsSoA = af::reorder(varsCopiedFromVec, 1, 2, 3, 0);
  for (int var=0; var<numVars; var++)
  {
    vars[var] = varsSoA(span, span, span, var);
  }
}

//void grid::setVecWithVars(array vars[])
//{
//
//}

void grid::setVarsWithVec(Vec vec)
{
  array varsCopiedFromVec;
  if (numGhost > 0)
  {
    DMGlobalToLocalBegin(dmGhost, vec, INSERT_VALUES, localVec);
    DMGlobalToLocalEnd(dmGhost, vec, INSERT_VALUES, localVec);

    double *pointerToLocalVec;
    VecGetArray(localVec, &pointerToLocalVec);

    varsCopiedFromVec = 
      array(numVars,
            gridParams::N1Local + 2*numGhost, 
            gridParams::N2Local + 2*numGhost,
            gridParams::N3Local + 2*numGhost,
            pointerToLocalVec
           );

    VecRestoreArray(globalVec, &pointerToLocalVec);
  }
  else
  {
    double *pointerToGlobalVec;
    VecGetArray(vec, &pointerToGlobalVec);

    varsCopiedFromVec = 
      array(numVars,
            gridParams::N1Local, gridParams::N2Local, gridParams::N3Local,
            pointerToGlobalVec
           );

    VecRestoreArray(vec, &pointerToGlobalVec);
  }

  array varsSoA = af::reorder(varsCopiedFromVec, 1, 2, 3, 0);
  for (int var=0; var<numVars; var++)
  {
    vars[var] = varsSoA(span, span, span, var);
  }

}

grid::~grid()
{
  delete [] vars;
  VecDestroy(&globalVec);
  VecDestroy(&localVec);

  DMDestroy(&dm);
  DMDestroy(&dmGhost);
}

