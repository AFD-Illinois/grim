#include "grim.h"

Vec conservedPetscVecOld;

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  
//  struct timeStepper ts;
//  timeStepperInit(&ts);
//  
//  timeStep(&ts);
//
//  timeStepperDestroy(&ts);


  DM dmdaWithoutGhostZones;

  Vec primPetscVec;
  Vec primPetscVecOld;
  Vec residualPetscVec;

  SNES snes;

  DMDACreate2d(PETSC_COMM_WORLD, 
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_BOX,
               N1, N2,
               PETSC_DECIDE, PETSC_DECIDE,
               DOF, 0, PETSC_NULL, PETSC_NULL, &dmdaWithoutGhostZones);

  SNESCreate(PETSC_COMM_WORLD, &snes);

  SNESSetDM(snes, dmdaWithoutGhostZones);

  DMCreateGlobalVector(dmdaWithoutGhostZones, &primPetscVec);
  DMCreateGlobalVector(dmdaWithoutGhostZones, &primPetscVecOld);
  DMCreateGlobalVector(dmdaWithoutGhostZones, &conservedPetscVecOld);
  DMCreateGlobalVector(dmdaWithoutGhostZones, &residualPetscVec);

  SNESSetFunction(snes, residualPetscVec, computeResidual, NULL);
  SNESSetFromOptions(snes);

  ARRAY(primGlobal);
  ARRAY(primOldGlobal);
  ARRAY(conservedOldGlobal);
  DMDAVecGetArrayDOF(dmdaWithoutGhostZones, 
                     primPetscVec, &primGlobal);
  DMDAVecGetArrayDOF(dmdaWithoutGhostZones, 
                     primPetscVecOld, &primOldGlobal);
  DMDAVecGetArrayDOF(dmdaWithoutGhostZones, conservedPetscVecOld, 
                     &conservedOldGlobal);

  int X1Start, X2Start, X3Start;
  int X1Size, X2Size, X3Size;

  DMDAGetCorners(dmdaWithoutGhostZones,
                 &X1Start, &X2Start, &X3Start,
                 &X1Size, &X2Size, &X3Size);

  PetscRandom random;
  PetscRandomCreate(PETSC_COMM_SELF, &random);
  PetscRandomSetType(random, PETSCRAND48);

  for (int j=X2Start; j<X2Start+X2Size; j++)
  {
    for (int i=X1Start; i<X1Start+X1Size; i++)
    {
      REAL XCoords[NDIM];

      REAL randomNum; 
      for (int var=0; var<DOF; var++)
      {
        PetscRandomGetValueReal(random, &randomNum);

        primOldGlobal[j][i][var] = randomNum;
      }

      for (int var=0; var<DOF; var++)
      {
        PetscRandomGetValueReal(random, &randomNum);

        primGlobal[j][i][var] = primOldGlobal[j][i][var] + 10.*randomNum;
      }

      struct geometry geom; setGeometry(XCoords, &geom);
      struct fluidElement elem;
      setFluidElement(primOldGlobal[j][i], &geom, &elem);

      REAL conservedVarsOld[DOF];
      computeFluxes(&elem, &geom, 0, conservedVarsOld);

      for (int var=0; var<DOF; var++)
      {
        conservedOldGlobal[j][i][var] = conservedVarsOld[var];
      }

    }
  }

  DMDAVecRestoreArrayDOF(dmdaWithoutGhostZones, 
                         primPetscVec, &primGlobal);
  DMDAVecRestoreArrayDOF(dmdaWithoutGhostZones, 
                         primPetscVecOld, &primOldGlobal);
  DMDAVecRestoreArrayDOF(dmdaWithoutGhostZones, conservedPetscVecOld, 
                         &conservedOldGlobal);


  SNESSolve(snes, NULL, primPetscVec);

  PetscFinalize();  
  return(0);
}


PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residualPetscVec,
                               void *ptr)
{
  int X1Start, X2Start, X3Start;
  int X1Size, X2Size, X3Size;


  DM dmdaWithoutGhostZones;
  SNESGetDM(snes, &dmdaWithoutGhostZones);

  DMDAGetCorners(dmdaWithoutGhostZones,
                 &X1Start, &X2Start, &X3Start,
                 &X1Size, &X2Size, &X3Size);

  ARRAY(primGlobal);
  ARRAY(conservedOldGlobal);
  ARRAY(residualGlobal);
  DMDAVecGetArrayDOF(dmdaWithoutGhostZones, primPetscVec, &primGlobal);
  DMDAVecGetArrayDOF(dmdaWithoutGhostZones, conservedPetscVecOld, 
                     &conservedOldGlobal);
  DMDAVecGetArrayDOF(dmdaWithoutGhostZones, residualPetscVec,
                     &residualGlobal);


  for (int j=X2Start; j<X2Start+X2Size; j++)
  {
    for (int i=X1Start; i<X1Start+X1Size; i++)
    {
      REAL XCoords[NDIM];

      struct geometry geom; setGeometry(XCoords, &geom);
      struct fluidElement elem;
      setFluidElement(primGlobal[j][i], &geom, &elem);

      REAL conservedVars[DOF];
      computeFluxes(&elem, &geom, 0, conservedVars);

      for (int var=0; var<DOF; var++)
      {
        residualGlobal[j][i][var] =   conservedVars[var] 
                                    - conservedOldGlobal[j][i][var];

      }
   
    }

  }

  DMDAVecRestoreArrayDOF(dmdaWithoutGhostZones, primPetscVec,
                         &primGlobal);
  DMDAVecRestoreArrayDOF(dmdaWithoutGhostZones, conservedPetscVecOld, 
                         &conservedOldGlobal);
  DMDAVecRestoreArrayDOF(dmdaWithoutGhostZones, residualPetscVec,
                         &residualGlobal);

  return(0);
}
