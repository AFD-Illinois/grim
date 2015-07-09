#include "timestepper.h"

/* Initialize a {\tt timeStepper} struct
 *
 * 1) 
 *
 * @param Input: ts, an uninitializied {\tt timeStepper} struct
 * @param Output: ts, an initialized {\tt timeStepper} struct
 */
void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1])
{
  /* If explicit or imex time stepping, then use a DM with no ghost zones. Graph
   * coloring with ensure that the number of evaluations to construct the
   * jacobian will be equal to DOF^2 */
  #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
    initGridData(N1, N2, N3, DOF, 0, &ts->primNPlusOne);
  #elif (TIME_STEPPING==IMPLICIT)
    initGridData(N1, N2, N3, DOF, NG, &ts->primNPlusOne);
  #endif

  initGridData(N1, N2, N3, DOF, NG, &ts->primNPlusHalf);
  initGridData(N1, N2, N3, DOF, NG, &ts->primN);

  initGridData(N1/TILE_SIZE_X1, N2, N3, DOF*TILE_SIZE_X1, 0, 
               &ts->conservedVarsN
               );
  initGridData(N1/TILE_SIZE_X1, N2, N3, DOF*TILE_SIZE_X1, 0, 
               &ts->divFluxes
               );
  initGridData(N1/TILE_SIZE_X1, N2, N3, DOF*TILE_SIZE_X1, 0, 
               &ts->sources
               );

  initGridData(N1, N2, N3, DOF, 0, &ts->residual);

  initGridData(N1/TILE_SIZE_X1, N2, N3, COMPUTE_DIM*TILE_SIZE_X1, 0, 
               &ts->dtGrid
              );
  initGridData(N1/TILE_SIZE_X1, N2, N3, 64*TILE_SIZE_X1, 0, 
               &ts->connection
              );

  initGridData(N1/TILE_SIZE_X1, N2, N3, 42*TILE_SIZE_X1, 0,
               &ts->geomCenter
              );
  initGridData(N1/TILE_SIZE_X1, N2, N3, 42*TILE_SIZE_X1, 0,
               &ts->geomFaceX1
              );
  initGridData(N1/TILE_SIZE_X1, N2, N3, 42*TILE_SIZE_X1, 0,
               &ts->geomFaceX1PlusOne
              );
  #if (COMPUTE_DIM == 2 || COMPUTE_DIM==3)
    initGridData(N1/TILE_SIZE_X1, N2, N3, 42*TILE_SIZE_X1, 0,
                 &ts->geomFaceX2
                );
    initGridData(N1/TILE_SIZE_X1, N2, N3, 42*TILE_SIZE_X1, 0,
                 &ts->geomFaceX2PlusOne
                );
  #endif
  #if (COMPUTE_DIM == 3)
    initGridData(N1/TILE_SIZE_X1, N2, N3, 42*TILE_SIZE_X1, 0,
                 &ts->geomFaceX3
                );
    initGridData(N1/TILE_SIZE_X1, N2, N3, 42*TILE_SIZE_X1, 0,
                 &ts->geomFaceX3PlusOne
                );
  #endif

  SNESCreate(PETSC_COMM_WORLD, &ts->snes);
  SNESSetDM(ts->snes, ts->primNPlusOne.dm);
  SNESSetFunction(ts->snes, ts->residual.vec, computeResidual, ts);
  SNESSetFromOptions(ts->snes);

  ts->dt = DT;
  ts->t = START_TIME;
  ts->tDump = START_TIME;

  ts->timeStepCounter = 0;
  ts->dumpCounter     = 0;

  ts->isZerothIterationOfSNES = 0;

  ts->computeDivOfFluxAtN         = 0;
  ts->computeDivOfFluxAtNPlusHalf = 0;
  ts->computeSourcesAtN           = 0;
  ts->computeSourcesAtNPlusHalf   = 0;

  /* Initialize problem dependent data */
  PetscMalloc1(1, &ts->problemSpecificData);

  if (ts->primN.iLocalSize % TILE_SIZE_X1 != 0)
  {
    PetscPrintf(PETSC_COMM_WORLD,
                "TILE_SIZE_X1 = %d does not divide iLocalSize = %d exactly\n",
                TILE_SIZE_X1, ts->primN.iLocalSize
               );
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  #if (USE_OPENMP)
  if (ts->primN.iLocalSize / TILE_SIZE_X1 % omp_get_max_threads() != 0)
  {
    PetscPrintf(PETSC_COMM_WORLD,
                "Tiles per OpenMP thread in X1 = %f not a perfect integer\n",
                  (REAL)ts->primN.iLocalSize \
                / (REAL) TILE_SIZE_X1 \
                / (REAL) omp_get_max_threads()
               );
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  #endif
  #if (COMPUTE_DIM==2 || COMPUTE_DIM==3)
    if (ts->primN.jLocalSize % TILE_SIZE_X2 != 0)
    {
      PetscPrintf(PETSC_COMM_WORLD,
                  "TILE_SIZE_X2 = %d does not divide jLocalSize = %d exactly\n",
                  TILE_SIZE_X2, ts->primN.jLocalSize
                 );
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    #if (USE_OPENMP)
    if (ts->primN.jLocalSize / TILE_SIZE_X2 % omp_get_max_threads() != 0)
    {
      PetscPrintf(PETSC_COMM_WORLD,
                  "Tiles per OpenMP thread in X2 = %f not a perfect integer\n",
                    (REAL)ts->primN.jLocalSize \
                  / (REAL) TILE_SIZE_X2 \
                  / (REAL) omp_get_max_threads()
                 );
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    #endif
  #endif

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD,
              "           Memory allocation complete\n\n");

  int numProcs;
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
  PetscPrintf(PETSC_COMM_WORLD,
              " Number of MPI procs being used       = %d\n",
              numProcs
             );
  #if (USE_OPENMP)
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of OpenMP threads being used  = %d\n",
                omp_get_max_threads()
               );
  #endif
  #if (COMPUTE_DIM==1)
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid size                            = %d\n",
                N1
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each MPI process      = %d\n",
                ts->primN.iLocalSize
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each tile             = %d\n",
                TILE_SIZE_X1
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each MPI process  = %d\n",
                ts->primN.iLocalSize/TILE_SIZE_X1
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each OpenMP thread= %d\n",
                ts->primN.iLocalSize/TILE_SIZE_X1/omp_get_max_threads()
               );
  #elif (COMPUTE_DIM==2)
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid size                            = %d x %d\n",
                N1, N2
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each MPI process      = %d x %d\n",
                ts->primN.iLocalSize, ts->primN.jLocalSize
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each tile             = %d x %d\n",
                TILE_SIZE_X1, TILE_SIZE_X2
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each MPI process  = %d x %d\n",
                ts->primN.iLocalSize/TILE_SIZE_X1,
                ts->primN.jLocalSize/TILE_SIZE_X2
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each OpenMP thread= %d x %d\n",
                ts->primN.iLocalSize/TILE_SIZE_X1/omp_get_max_threads(),
                ts->primN.jLocalSize/TILE_SIZE_X2/omp_get_max_threads()
               );
  #elif (COMPUTE_DIM==3)
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid size                            = %d x %d x %d\n",
                N1, N2, N3
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each MPI process      = %d x %d x %d\n",
                ts->primN.iLocalSize, ts->primN.jLocalSize, 
                ts->primN.kLocalSize
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Grid points in each tile             = %d x %d x %d\n",
                TILE_SIZE_X1, TILE_SIZE_X2, 2*NG + 1
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each MPI process  = %d x %d\n",
                ts->primN.iLocalSize/TILE_SIZE_X1,
                ts->primN.jLocalSize/TILE_SIZE_X2
               );
    PetscPrintf(PETSC_COMM_WORLD,
                " Number of tiles in each OpenMP thread= %d x %d\n",
                ts->primN.iLocalSize/TILE_SIZE_X1/omp_get_max_threads(),
                ts->primN.jLocalSize/TILE_SIZE_X2/omp_get_max_threads()
               );
  #endif

  PetscPrintf(PETSC_COMM_WORLD,
              "#################################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  /* Precompute the Chritoffel symbols gammaUpDownDown */
//  PetscPrintf(PETSC_COMM_WORLD, "Computing Christoffel symbols...");
//  setChristoffelSymbols(ts);
//  PetscPrintf(PETSC_COMM_WORLD, "done\n");

  #if (RESTART)
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if (rank==0)
    {
      if (access(RESTART_FILE, F_OK) != -1)
      {
        /* File exists */
        PetscPrintf(PETSC_COMM_WORLD, "\nFound restart file: %s\n\n", RESTART_FILE); 
      }
      else
      {
        /* File does not exist */
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        PetscPrintf(PETSC_COMM_WORLD, "Restart file %s does not exist\n",
                    RESTART_FILE
                   );
        MPI_Abort(PETSC_COMM_WORLD, 1);
      }
    }

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, "restartfile.h5",
                        FILE_MODE_READ, &viewer);
    PetscObjectSetName((PetscObject) ts->primPetscVecOld, "primVars");
    VecLoad(ts->primPetscVecOld, viewer);
    PetscViewerDestroy(&viewer);

  #else

    /* Set initialConditions from problem */
    initialConditions(ts);

  #endif /* RESTART option */

  /* Output the initial conditions */
//  VecCopy(ts->primPetscVecOld, ts->primPetscVec);
//  diagnostics(ts);

}

//#if (CONDUCTION)
//void initConductionDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
//{
//
//  #if (COMPUTE_DIM==1)
//    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM, 0, NULL,
//                 &ts->gradTDM);
//
//    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM*NDIM, 0, NULL,
//                 &ts->graduConDM);
//
//    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM, 0, NULL,
//                 &ts->graduConHigherOrderTermsDM);
//
//  #elif (COMPUTE_DIM==2)
//    DMDACreate2d(PETSC_COMM_WORLD, 
//                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
//                 DMDA_STENCIL_BOX,
//                 N1, N2,
//                 PETSC_DECIDE, PETSC_DECIDE,
//                 COMPUTE_DIM, 0, PETSC_NULL, PETSC_NULL,
//                 &ts->gradTDM);
//
//    DMDACreate2d(PETSC_COMM_WORLD, 
//                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
//                 DMDA_STENCIL_BOX,
//                 N1, N2,
//                 PETSC_DECIDE, PETSC_DECIDE,
//                 COMPUTE_DIM*NDIM, 0, PETSC_NULL, PETSC_NULL,
//                 &ts->graduConDM);
//
//    DMDACreate2d(PETSC_COMM_WORLD, 
//                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
//                 DMDA_STENCIL_BOX,
//                 N1, N2,
//                 PETSC_DECIDE, PETSC_DECIDE,
//                 COMPUTE_DIM, 0, PETSC_NULL, PETSC_NULL,
//                 &ts->graduConHigherOrderTermsDM);
//
//  #endif
//
//  DMCreateGlobalVector(ts->gradTDM, &ts->gradTPetscVec);
//  DMCreateGlobalVector(ts->graduConDM, &ts->graduConPetscVec);
//  DMCreateGlobalVector(ts->graduConHigherOrderTermsDM, 
//                       &ts->graduConHigherOrderTerm1PetscVec);
//  DMCreateGlobalVector(ts->graduConHigherOrderTermsDM, 
//                       &ts->graduConHigherOrderTerm2PetscVec);
//
//  VecSet(ts->gradTPetscVec, 0.);
//  VecSet(ts->graduConPetscVec, 0.);
//  VecSet(ts->graduConHigherOrderTerm1PetscVec, 0.);
//  VecSet(ts->graduConHigherOrderTerm2PetscVec, 0.);
//}
//
//void destroyConductionDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
//{
//  VecDestroy(&ts->gradTPetscVec);
//  VecDestroy(&ts->graduConPetscVec);
//  VecDestroy(&ts->graduConHigherOrderTerm1PetscVec);
//  VecDestroy(&ts->graduConHigherOrderTerm2PetscVec);
//
//  DMDestroy(&ts->gradTDM);
//  DMDestroy(&ts->graduConDM);
//  DMDestroy(&ts->graduConHigherOrderTermsDM);
//}
//#endif
//
//
//#if (VISCOSITY)
//void initViscosityDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
//{
//
//  #if (COMPUTE_DIM==1)
//    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM*NDIM, 0, NULL,
//                 &ts->graduConVisDM);
//    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, COMPUTE_DIM, 0, NULL,
//                 &ts->graduConHigherOrderTermsVisDM);
//
//  #elif (COMPUTE_DIM==2)
//    DMDACreate2d(PETSC_COMM_WORLD, 
//                   DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
//                   DMDA_STENCIL_BOX,
//                   N1, N2,
//                   PETSC_DECIDE, PETSC_DECIDE,
//                   COMPUTE_DIM*NDIM, 0, PETSC_NULL, PETSC_NULL,
//                   &ts->graduConVisDM);
//    DMDACreate2d(PETSC_COMM_WORLD, 
//                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
//                 DMDA_STENCIL_BOX,
//                 N1, N2,
//                 PETSC_DECIDE, PETSC_DECIDE,
//                 COMPUTE_DIM, 0, PETSC_NULL, PETSC_NULL,
//                 &ts->graduConHigherOrderTermsVisDM);
//
//  #endif
//
//  DMCreateGlobalVector(ts->graduConVisDM, &ts->graduConVisPetscVec);
//  DMCreateGlobalVector(ts->graduConHigherOrderTermsVisDM, 
//                       &ts->graduConHigherOrderTerm1VisPetscVec);
//  DMCreateGlobalVector(ts->graduConHigherOrderTermsVisDM, 
//                       &ts->graduConHigherOrderTerm2VisPetscVec);
//
//  VecSet(ts->graduConVisPetscVec, 0.);
//  VecSet(ts->graduConHigherOrderTerm1VisPetscVec, 0.);
//  VecSet(ts->graduConHigherOrderTerm2VisPetscVec, 0.);
//}
//
//void destroyViscosityDataStructures(struct timeStepper ts[ARRAY_ARGS 1])
//{
//  VecDestroy(&ts->graduConVisPetscVec);
//  VecDestroy(&ts->graduConHigherOrderTerm1VisPetscVec);
//  VecDestroy(&ts->graduConHigherOrderTerm2VisPetscVec);
//
//  DMDestroy(&ts->graduConVisDM);
//  DMDestroy(&ts->graduConHigherOrderTermsVisDM);
//}
//#endif

void setChristoffelSymbols(struct timeStepper ts[ARRAY_ARGS 1])
{
  setPointerToVec(&ts->connection);

  LOOP_OVER_TILES(&ts->connection)
  {
    struct gridTile tile;
    setTile(iTile, jTile, kGlobal, &ts->connection, &tile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2, 0, TILE_SIZE_X3)
    {
      
      struct gridZone zone;
      setZone(iInTile, jInTile, kInTile, &tile, &zone);

      REAL XCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      struct geometry geom; setGeometry(XCoords, &geom);

      /* Now compute connectionGlobal with 
       * Index Up   - eta
       * Index down - mu
       * Index down - nu */
      for (int eta=0; eta<NDIM; eta++)
      {
        for (int mu=0; mu<NDIM; mu++)
        {
          for (int nu=0; nu<NDIM; nu++)
          {
            for (int alpha=0; alpha<NDIM; alpha++)
            {
              INDEX_GRID(&ts->connection, &zone, GAMMA_UP_DOWN_DOWN(eta, mu, nu))
                +=
                  geom.gCon[eta][alpha]
                * gammaDownDownDown(alpha, mu, nu, XCoords);
            }
          }
        }
      }

    }
  }

  restorePointerToVec(&ts->connection);
}

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1])
{
  destroyGridData(&ts->primNPlusOne);
  destroyGridData(&ts->primNPlusHalf);
  destroyGridData(&ts->primN);

  destroyGridData(&ts->conservedVarsN);
  destroyGridData(&ts->divFluxes);
  destroyGridData(&ts->sources);

  destroyGridData(&ts->residual);
  destroyGridData(&ts->dtGrid);
  destroyGridData(&ts->connection);

  SNESDestroy(&ts->snes);

  PetscFree(ts->problemSpecificData);

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "# Memory deallocation complete #\n");
  PetscPrintf(PETSC_COMM_WORLD, "################################\n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");
}
