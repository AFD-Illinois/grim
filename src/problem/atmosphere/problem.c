#include "../problem.h"

#if (CONDUCTION)
REAL kappaProblem;
REAL tauProblem;

void setConductionParameters(struct fluidElement elem[ARRAY_ARGS 1])
{
  elem->kappa = kappaProblem;
  elem->tau   = tauProblem;
}
#endif

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  REAL randNum;
  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);

  REAL rho1D[N1+2*NG];
  REAL u1D[N1+2*NG];
  REAL uUpr1D[N1+2*NG];
  REAL phi1D[N1+2*NG];
  REAL rCoords1D[N1+2*NG];

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); /* get current proc id */

  if (rank==0)
  {
    FILE *rhoFile;
    FILE *uFile;
    FILE *uUprFile;
    FILE *phiFile;
    FILE *rCoordsFile;
  
    char *rhoLine     = NULL;
    char *uLine       = NULL;
    char *uUprLine    = NULL;
    char *phiLine     = NULL; 
    char *rCoordsLine = NULL;

    size_t rhoLen     = 0; ssize_t rhoRead;
    size_t uLen       = 0; ssize_t uRead;
    size_t uUprLen    = 0; ssize_t uUprRead;
    size_t phiLen     = 0; ssize_t phiRead;
    size_t rCoordsLen = 0; ssize_t rCoordsRead;

    rhoFile     = fopen(RHO_INPUT_FILE    , "r");
    uFile       = fopen(UU_INPUT_FILE     , "r");
    uUprFile    = fopen(UR_INPUT_FILE     , "r");
    phiFile     = fopen(PHI_INPUT_FILE    , "r");
    rCoordsFile = fopen(RCOORDS_INPUT_FILE, "r");

    if (   rhoFile      == NULL 
        || uFile        == NULL
        || uUprFile     == NULL
        || phiFile      == NULL
        || rCoordsFile  == NULL
       )
    {
      SETERRQ(PETSC_COMM_WORLD, 1, "Input data files not found!\n");
    }

    for (int i=-NG; i<N1+NG; i++)
    {
      rhoRead     = getline(&rhoLine    , &rhoLen     , rhoFile);
      uRead       = getline(&uLine      , &uLen       , uFile);
      uUprRead    = getline(&uUprLine   , &uUprLen    , uUprFile);
      phiRead     = getline(&phiLine    , &phiLen     , phiFile);
      rCoordsRead = getline(&rCoordsLine, &rCoordsLen , rCoordsFile);

      if (   rhoRead   == -1
          || uRead     == -1
          || uUprRead  == -1
          || phiRead   == -1
          || rCoordsRead == -1
         )
      {
        SETERRQ(PETSC_COMM_WORLD, 1, 
                "Found the input data files but failed to read them!\n");
      }

      rho1D[i+NG]     = atof(rhoLine);
      u1D[i+NG]       = atof(uLine);
      uUpr1D[i+NG]    = atof(uUprLine);
      phi1D[i+NG]     = atof(phiLine);
      rCoords1D[i+NG] = atof(rCoordsLine);
    }

    free(rhoLine);
    free(uLine);
    free(uUprLine);
    free(phiLine);
    free(rCoordsLine);

    fclose(rhoFile);
    fclose(uFile);
    fclose(uUprFile);
    fclose(phiFile);
    fclose(rCoordsFile);
  }

  /* Broadcast the data from rank 0 proc to all other procs */
  MPI_Bcast(&rho1D[0]     , N1+2*NG, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&u1D[0]       , N1+2*NG, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&uUpr1D[0]    , N1+2*NG, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&phi1D[0]     , N1+2*NG, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&rCoords1D[0] , N1+2*NG, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  /* Broadcast complete */

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM], xCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      XTox(XCoords, xCoords);

      REAL r = xCoords[1];
      if (fabs(r - rCoords1D[zone.i+NG]) > 1e-15)
      {
        SETERRQ2(PETSC_COMM_WORLD, 1, 
                 "Mismatch in rCoords! Check r coords in python script. r = %f, rCoords = %f\n",
                 r, rCoords1D[zone.i+NG]
                );
      }

      primTile[INDEX_TILE(&zone, RHO)] = rho1D[zone.i+NG];
      PetscRandomGetValue(randNumGen, &randNum);
      primTile[INDEX_TILE(&zone, UU)]  = 
        u1D[zone.i+NG]*(1. + PERTURBATIONS_AMPLITUDE*(randNum-0.5));

      #if (CONDUCTION)
        primTile[INDEX_TILE(&zone, PHI)] = phi1D[zone.i+NG];
        
        kappaProblem = 0.1;
        tauProblem   = 10.;
      #endif

      REAL uConBL[NDIM];
      uConBL[1] = uUpr1D[zone.i+NG];
      uConBL[2] = 0.;
      uConBL[3] = 0.;

      /* Formula to output vUpr in MKS from uUpr in BL from Ben Ryan */
      struct geometry geom;
      setGeometry(XCoords, &geom);
      REAL a = geom.gCov[1][1];
      REAL b = geom.gCon[0][1];
      REAL c = geom.gCon[0][0];
      REAL vUprMKS = 
        (  c*uConBL[1]/r 
         - sqrt(-a*b*b*b*b
                -a*b*b*c*uConBL[1]*uConBL[1]/(r*r) - b*b*c
               )
        )/(a*b*b + c);

      primTile[INDEX_TILE(&zone, U1)] = vUprMKS;
      primTile[INDEX_TILE(&zone, U2)] = 0.; 
      primTile[INDEX_TILE(&zone, U3)] = 0.; 

      primTile[INDEX_TILE(&zone, B1)] = 0.;
      primTile[INDEX_TILE(&zone, B2)] = 0.;
      primTile[INDEX_TILE(&zone, B3)] = 0.;

      /* Set up Dirichlet boundary data */

      if (zone.i < 0)
      {
        for (int var=0; var<DOF; var++)
        {
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][var]
            = primTile[INDEX_TILE(&zone, var)];
        }
      }

      if (zone.i > N1-1)
      {
        for (int var=0; var<DOF; var++)
        {
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][var]
            = primTile[INDEX_TILE(&zone, var)];
        }
      }
  
    }

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);
      
      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primOldGlobal, &zone, var) =
          primTile[INDEX_TILE(&zone, var)];
      }
    }

  }


  Vec bSqrPetscVecGlobal;
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &bSqrPetscVecGlobal);

  ARRAY(bSqrGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, bSqrPetscVecGlobal,
                     &bSqrGlobal);


  /* Now set the magnetic field using the magnetic vector potential */
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL AVectorTile[TILE_SIZE];

    LOOP_INSIDE_TILE(-1, TILE_SIZE_X1+1, -1, TILE_SIZE_X2+1)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM], xCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      XTox(XCoords, xCoords);

      REAL r     = xCoords[1];
      REAL theta = xCoords[2];

      REAL AVec = 0.;
      #if (MAGNETIC_FIELD_CONFIGURATION == VERTICAL)

        AVec = r*sin(theta);

      #endif

      AVectorTile[INDEX_TILE(&zone, 0)] = AVec;

    }

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      struct geometry geom; setGeometry(XCoords, &geom);
      REAL g = sqrt(-geom.gDet);

      INDEX_PETSC(primOldGlobal, &zone, B1) = 
        -(  AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile,   0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile+1, 0)]
          + AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile,   0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile+1, 0)]
         )/(2.*zone.dX2*g);

      INDEX_PETSC(primOldGlobal, &zone, B2) = 
         (  AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile,   0)]
          + AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile+1, 0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile,   0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile+1, 0)]
         )/(2.*zone.dX1*g);

      INDEX_PETSC(primOldGlobal, &zone, B3) = 0.; 

      #if (MAGNETIC_FIELD_CONFIGURATION == THETA)
        INDEX_PETSC(primOldGlobal, &zone, B1) = 0.;
        INDEX_PETSC(primOldGlobal, &zone, B2) = 1./g;
        INDEX_PETSC(primOldGlobal, &zone, B3) = 0.;
      #elif (MAGNETIC_FIELD_CONFIGURATION == RADIAL)
        INDEX_PETSC(primOldGlobal, &zone, B1) = 1./g;
        INDEX_PETSC(primOldGlobal, &zone, B2) = 0.;
        INDEX_PETSC(primOldGlobal, &zone, B3) = 0.;
      #endif

      struct fluidElement elem;
      setFluidElement(&INDEX_PETSC(primOldGlobal, &zone, 0), &geom, &elem);

      REAL bCov[NDIM];
      conToCov(elem.bCon, &geom, bCov);
      REAL bSqr = covDotCon(bCov, elem.bCon);

      INDEX_PETSC(bSqrGlobal, &zone, 0) = bSqr;
    } 

  }

  REAL bSqrMax, uMax;
  VecStrideMax(bSqrPetscVecGlobal, 0, NULL, &bSqrMax);
  VecStrideMax(ts->primPetscVecOld, UU, NULL, &uMax);

  REAL betaActual = (ADIABATIC_INDEX-1.)*uMax/(0.5*bSqrMax);
  REAL norm = sqrt(betaActual/PLASMA_BETA);

  VecStrideScale(ts->primPetscVecOld, B1, norm);
  VecStrideScale(ts->primPetscVecOld, B2, norm);
  VecStrideScale(ts->primPetscVecOld, B3, norm);

  VecDestroy(&bSqrPetscVecGlobal);

  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, bSqrPetscVecGlobal,
                         &bSqrGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      struct geometry geom; setGeometry(XCoords, &geom);
      REAL g = sqrt(-geom.gDet);

      if (zone.i < 0)
      {
        #if (MAGNETIC_FIELD_CONFIGURATION == RADIAL)
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][B1]
            = norm/g;
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][B2]
            = 0.;
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][B3]
            = 0.;
        #elif (MAGNETIC_FIELD_CONFIGURATION == THETA)
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][B1]
            = 0.;
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][B2]
            = norm/g;
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][B3]
            = 0.;
        #endif
      }

      if (zone.i > N1-1)
      {
        #if (MAGNETIC_FIELD_CONFIGURATION == RADIAL)
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][B1]
            = norm/g;
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][B2]
            = 0.;
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][B3]
            = 0.;
        #elif (MAGNETIC_FIELD_CONFIGURATION == THETA)
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][B1]
            = 0.;
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][B2]
            = norm/g;
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][B3]
            = 0.;
        #endif
      }
    } 
  }

  PetscPrintf(PETSC_COMM_WORLD, "done\n");
  /* Done with setting the initial conditions */
}

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
                                     const int X1Start, const int X2Start,
                                     const int X1Size, const int X2Size,
                                     REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    if (primTile[INDEX_TILE(&zone, RHO)] < RHO_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR_MIN;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < UU_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR_MIN;
    }
  }
}

void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    REAL XCoords[NDIM], xCoords[NDIM];
    getXCoords(&zone, CENTER, XCoords);
    XTox(XCoords, xCoords);

    REAL r = xCoords[1];

    REAL rhoFloor = RHO_FLOOR*pow(r, RHO_FLOOR_FALLOFF);
    REAL uFloor = UU_FLOOR*pow(r, UU_FLOOR_FALLOFF);

    if (rhoFloor < RHO_FLOOR_MIN)
    {
      rhoFloor = RHO_FLOOR_MIN;
    }

    if (uFloor < UU_FLOOR_MIN)
    {
      uFloor = UU_FLOOR_MIN;
    }

    REAL rho = primTile[INDEX_TILE(&zone, RHO)];
    REAL u = primTile[INDEX_TILE(&zone, UU)];

    if (rho < rhoFloor)
    {
      primTile[INDEX_TILE(&zone, RHO)] = rhoFloor;
    }

    if (u < uFloor)
    {
      primTile[INDEX_TILE(&zone, UU)] = uFloor;
    }

  }

}

void applyAdditionalProblemSpecificBCs
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL primTile[ARRAY_ARGS TILE_SIZE]
)
{

  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    #if (PHYSICAL_BOUNDARY_LEFT_EDGE==DIRICHLET)

      #if (MAGNETIC_FIELD_CONFIGURATION == VERTICAL)
        SETERRQ(PETSC_COMM_WORLD, 1,
                "Dirichet boundary conditions not implemented for \
                 vertical magnetic field configuration\n");
      #endif

      if (zone.i < 0)
      {
        for (int var=0; var<DOF; var++)
        {
          primTile[INDEX_TILE(&zone, var)] = 
            problemSpecificData->primVarsLeftEdge[zone.i+NG][var];
        }
      }
    #endif

    #if (PHYSICAL_BOUNDARY_RIGHT_EDGE==DIRICHLET)
      
      #if (MAGNETIC_FIELD_CONFIGURATION == VERTICAL)
        SETERRQ(PETSC_COMM_WORLD, 1,
                "Dirichet boundary conditions not implemented for \
                 vertical magnetic field configuration\n");
      #endif

      if (zone.i > N1-1)
      {
        for (int var=0; var<DOF; var++)
        {
          primTile[INDEX_TILE(&zone, var)] = 
            problemSpecificData->primVarsRightEdge[zone.i-N1][var];
        }
      }
    #endif

  }

}

void applyProblemSpecificFluxFilter
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
)
{

}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primHalfStepGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecHalfStep, &primHalfStepGlobal);

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] =
        INDEX_PETSC(primHalfStepGlobal, &zone, var);
      }
    }

    applyFloor(iTile, jTile,
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primHalfStepGlobal, &zone, var) =
        primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecHalfStep, &primHalfStepGlobal);
}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] =
        INDEX_PETSC(primOldGlobal, &zone, var);
      }
    }

    applyFloor(iTile, jTile,
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primOldGlobal, &zone, var) =
        primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);
}
