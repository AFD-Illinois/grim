#include "grid.hpp"

grid::grid(int numVars)
{
  this->numVars  = numVars;
  this->numGhost = params::numGhost; 

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
                   params::N1, numVars, numGhost, NULL,
                   &dm
                  );

      domainX1 = new af::seq(params::numGhost, af::end - params::numGhost - 1);
      domainX2 = new af::seq(span);
      domainX3 = new af::seq(span);
      break;
  
    case 2:

      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost,
                   PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      domainX1 = new af::seq(params::numGhost, af::end - params::numGhost - 1);
      domainX2 = new af::seq(params::numGhost, af::end - params::numGhost - 1);
      domainX3 = new af::seq(span);
      break;

    case 3:

      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2, params::N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost, 
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dm
                  );
      domainX1 = new af::seq(params::numGhost, af::end - params::numGhost - 1);
      domainX2 = new af::seq(params::numGhost, af::end - params::numGhost - 1);
      domainX3 = new af::seq(params::numGhost, af::end - params::numGhost - 1);
      break;
  }

  DMCreateGlobalVector(dm, &globalVec);
  DMCreateLocalVector(dm, &localVec);
  VecSet(globalVec, 0.);
  VecSet(localVec,  0.);

  DMDAGetCorners
    (dm, &iLocalStart, &jLocalStart, &kLocalStart,
         &N1Local,     &N2Local,     &N3Local
    );

  iLocalEnd = iLocalStart + N1Local;
  jLocalEnd = jLocalStart + N2Local;
  kLocalEnd = kLocalStart + N3Local;

  dX1 = (params::X1End - params::X1Start)/params::N1;
  dX2 = (params::X2End - params::X2Start)/params::N2;
  dX3 = (params::X3End - params::X3Start)/params::N3;

  vars = new array[numVars];

  /* Initialize vars to localVec */
  copyLocalVecToVars();
}

void grid::copyLocalVecToVars()
{
  array varsCopiedFromLocalVec;

  double *pointerToLocalVec;
  VecGetArray(localVec, &pointerToLocalVec);

  switch (params::dim)
  {
    case 1:
      varsCopiedFromLocalVec = 
        array(numVars,
              N1Local + 2*numGhost, 
              N2Local,
              N3Local,
              pointerToLocalVec
             );

      break;

    case 2:
      varsCopiedFromLocalVec = 
        array(numVars,
              N1Local + 2*numGhost, 
              N2Local + 2*numGhost,
              N3Local,
              pointerToLocalVec
             );

      break;

    case 3:
      varsCopiedFromLocalVec = 
        array(numVars,
              N1Local + 2*numGhost, 
              N2Local + 2*numGhost,
              N3Local + 2*numGhost,
              pointerToLocalVec
             );

      break;
    }

  VecRestoreArray(localVec, &pointerToLocalVec);

  varsSoA = af::reorder(varsCopiedFromLocalVec, 1, 2, 3, 0);
  for (int var=0; var<numVars; var++)
  {
    vars[var] = varsSoA(span, span, span, var);
  }
}

void grid::communicate()
{
  /* Get data into Array of Structs format as needed by Petsc */
  array varsAoS = af::reorder(varsSoA, 3, 0, 1, 2);

  /* Copy part of varsAoS contained in [0, N1) x [0, N2) x [0, N3) to global vec */
  double *pointerToLocalVec;
  double *pointerToGlobalVec;
  VecGetArray(localVec, &pointerToLocalVec);
  VecGetArray(globalVec, &pointerToGlobalVec);

  varsAoS.host(pointerToLocalVec);

  for (int k=0; k<N3Local; k++)
  {
    for (int j=0; j<N2Local; j++)
    {
      for (int i=0; i<N1Local; i++)
      {
        for (int var=0; var<numVars; var++)
        {
          int index = var + numVars*(i + N1Local*(j + N2Local*(k) ) );
          pointerToGlobalVec[index] = pointerToLocalVec[index];
        }
      }
    }
  }

  VecRestoreArray(localVec, &pointerToLocalVec);
  VecRestoreArray(globalVec, &pointerToGlobalVec);

  DMGlobalToLocalBegin(dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(dm, globalVec, INSERT_VALUES, localVec);

  /* Now copy data from localVec to vars */
  copyLocalVecToVars();
}



grid::~grid()
{
  delete [] vars;
  VecDestroy(&globalVec);
  VecDestroy(&localVec);

  DMDestroy(&dm);
}

