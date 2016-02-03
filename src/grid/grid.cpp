#include "grid.hpp"

grid::grid(int N1, int N2, int N3, 
           int numGhost, int dim, int numVars,
           DMBoundaryType boundaryLeft,
           DMBoundaryType boundaryRight,
           DMBoundaryType boundaryTop,
           DMBoundaryType boundaryBottom,
           DMBoundaryType boundaryFront,
           DMBoundaryType boundaryBack
          )
{
  this->numVars  = numVars;
  this->numGhost = numGhost;
  this->N1 = N1;
  this->N2 = N2;
  this->N3 = N3;
  this->boundaryLeft   = boundaryLeft;
  this->boundaryRight  = boundaryRight;
  this->boundaryTop    = boundaryTop;
  this->boundaryBottom = boundaryBottom;
  this->boundaryFront  = boundaryFront;
  this->boundaryBack   = boundaryBack;

  hasHostPtrBeenAllocated = 0;

  /* Implementations for MIRROR, OUTFLOW in boundary.cpp and DIRICHLET in
   * problem.cpp */
  boundaryLeft  = DM_BOUNDARY_GHOSTED; boundaryRight  = DM_BOUNDARY_GHOSTED;
  boundaryTop   = DM_BOUNDARY_GHOSTED; boundaryBottom = DM_BOUNDARY_GHOSTED;
  boundaryFront = DM_BOUNDARY_GHOSTED; boundaryBack   = DM_BOUNDARY_GHOSTED;

  if (   boundaryLeft  == boundaries::PERIODIC 
      || boundaryRight == boundaries::PERIODIC
     )
  {
    boundaryLeft  = DM_BOUNDARY_PERIODIC;
    boundaryRight = DM_BOUNDARY_PERIODIC;
  }

  if (   boundaryTop    == boundaries::PERIODIC 
      || boundaryBottom == boundaries::PERIODIC
     )
  {
    boundaryTop    = DM_BOUNDARY_PERIODIC;
    boundaryBottom = DM_BOUNDARY_PERIODIC;
  }

  if (   boundaryFront == boundaries::PERIODIC 
      || boundaryBack  == boundaries::PERIODIC
     )
  {
    boundaryFront = DM_BOUNDARY_PERIODIC;
    boundaryBack  = DM_BOUNDARY_PERIODIC;
  }

  switch (params::dim)
  {
    case 1:
      N1Total = N1 + 2*numGhost;
      N2Total = 1;
      N3Total = 1;

      numGhostX1 = numGhost;
      numGhostX2 = 0;
      numGhostX3 = 0;

      domainX1 = new af::seq(numGhost, af::end - numGhost - 1);
      domainX2 = new af::seq(span);
      domainX3 = new af::seq(span);

      DMDACreate1d(PETSC_COMM_WORLD, boundaryLeft, 
                   N1, numVars, numGhostX1, NULL,
                   &dm
                  );
      break;
  
    case 2:
      N1Total = N1 + 2*numGhost;
      N2Total = N2 + 2*numGhost;
      N3Total = 1;

      numGhostX1 = numGhost;
      numGhostX2 = numGhost;
      numGhostX3 = 0;

      domainX1 = new af::seq(numGhost, af::end - numGhost - 1);
      domainX2 = new af::seq(numGhost, af::end - numGhost - 1);
      domainX3 = new af::seq(span);

      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   N1, N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhostX1,
                   PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      break;

    case 3:
      N1Total = N1 + 2*numGhost;
      N2Total = N2 + 2*numGhost;
      N3Total = N3 + 2*numGhost;

      numGhostX1 = numGhost;
      numGhostX2 = numGhost;
      numGhostX3 = numGhost;

      domainX1 = new af::seq(numGhost, af::end - numGhost - 1);
      domainX2 = new af::seq(numGhost, af::end - numGhost - 1);
      domainX3 = new af::seq(numGhost, af::end - numGhost - 1);

      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   N1, N2, N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhostX1,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dm
                  );

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

  varsCopiedFromLocalVec = 
    array(numVars,
          N1Total, N2Total, N3Total, 
          pointerToLocalVec
         );

  VecRestoreArray(localVec, &pointerToLocalVec);

  varsSoA = af::reorder(varsCopiedFromLocalVec, 1, 2, 3, 0);
  for (int var=0; var<numVars; var++)
  {
    vars[var] = varsSoA(span, span, span, var);
  }
}

void grid::copyVarsToHostPtr()
{
  for (int var=0; var < numVars; var++)
  {
            /*i,    j,    k*/
    varsSoA(span, span, span, var) = vars[var];
  }

  /* If already called once, free the existing memory */
  if (hasHostPtrBeenAllocated)
  {
    delete hostPtr;
  }

  /* Allocate space for hostPtr and then copy data from device to host */
  hostPtr =  varsSoA.host<double>();
  hasHostPtrBeenAllocated = 1;
}

void grid::copyHostPtrToVars(const double *hostPtr)
{
  varsSoA = array(N1Total, N2Total, N3Total, numVars,
                  hostPtr
                 );

  for (int var=0; var < numVars; var++)
  {
    vars[var] = varsSoA(span, span, span, var);
  }
}

void grid::communicate()
{
  /* Get data into Array of Structs format as needed by Petsc */
  for (int var=0; var < numVars; var++)
  {
    varsSoA(span, span, span, var) = vars[var];
  }

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
	        /* Note on indices: Petsc uses non-ghosted vectors for the global vector,
	        * but ghosted indices for the local vector! */
          const int globalindex = var + numVars*(i + N1Local*(j + N2Local*(k) ) );
          const int localindex  = var + numVars*(i + numGhostX1 
                                                   + N1Total*(j + numGhostX2
                                                              + N2Total*(k + numGhostX3)
                                                             ) 
                                                );

          pointerToGlobalVec[globalindex] = pointerToLocalVec[localindex];
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
  if (hasHostPtrBeenAllocated)
  {
    delete hostPtr;
  }
  delete [] vars;
  VecDestroy(&globalVec);
  VecDestroy(&localVec);

  DMDestroy(&dm);

  delete domainX1, domainX2, domainX3;
}
