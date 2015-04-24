#include "timestepper.h"

void diagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  if (ts->timeStepCounter==0)
  {
    /* Dump geometry data */

    PetscPrintf(PETSC_COMM_WORLD, "Dumping the metric and the Christoffel symbols...");

    DM metricDMDA;
    Vec gCovPetscVec, gConPetscVec;

#if (COMPUTE_DIM==1)
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, 16, 0, NULL,
                 &metricDMDA);
#elif (COMPUTE_DIM==2)
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX, N1, N2, PETSC_DECIDE, PETSC_DECIDE,
                 16, 0, PETSC_NULL, PETSC_NULL, &metricDMDA);
#endif /* Choose dimension */

    DMCreateGlobalVector(metricDMDA, &gCovPetscVec);
    DMCreateGlobalVector(metricDMDA, &gConPetscVec);

    ARRAY(gCovGlobal);
    ARRAY(gConGlobal);
    DMDAVecGetArrayDOF(metricDMDA, gCovPetscVec, &gCovGlobal);
    DMDAVecGetArrayDOF(metricDMDA, gConPetscVec, &gConGlobal);

    LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
    {
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

        struct geometry geom;
        setGeometry(XCoords, &geom);

        for (int mu=0; mu<NDIM; mu++)
        {
          for (int nu=0; nu<NDIM; nu++)
          {
            #if (COMPUTE_DIM==1)
              gCovGlobal[zone.i][nu + NDIM*mu] = geom.gCov[mu][nu];
              gConGlobal[zone.i][nu + NDIM*mu] = geom.gCon[mu][nu];
            #elif (COMPUTE_DIM==2)
              gCovGlobal[zone.j][zone.i][nu + NDIM*mu] = geom.gCov[mu][nu];
              gConGlobal[zone.j][zone.i][nu + NDIM*mu] = geom.gCon[mu][nu];
            #endif
          }
        }

      }
    }

    DMDAVecRestoreArrayDOF(metricDMDA, gCovPetscVec, &gCovGlobal);
    DMDAVecRestoreArrayDOF(metricDMDA, gConPetscVec, &gConGlobal);

    char parametersFilename[50];
    sprintf(parametersFilename, "parameters.h5");

    PetscViewer parametersViewer;
    
    PetscViewerHDF5Open(PETSC_COMM_WORLD, parametersFilename,
                        FILE_MODE_WRITE, &parametersViewer);

    PetscObjectSetName((PetscObject) gCovPetscVec, "gCov");
    PetscObjectSetName((PetscObject) gConPetscVec, "gCon");
    PetscObjectSetName((PetscObject) ts->connectionPetscVec,
                       "gammaUpdowndown");

    VecView(gCovPetscVec, parametersViewer);
    VecView(gConPetscVec, parametersViewer);
    VecView(ts->connectionPetscVec, parametersViewer);

    /* Output parameters */
    /* Available data types:
       http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscDataType.html
         PETSC_INT
         PETSC_DOUBLE
         PETSC_COMPLEX
         PETSC_LONG
         PETSC_SHORT
         PETSC_FLOAT
         PETSC_CHAR
         PETSC_BIT_LOGICAL
         PETSC_ENUM
         PETSC_BOOL
         PETSC___FLOAT128
         PETSC_OBJECT
         PETSC_FUNCTION
         PETSC_STRING
     */

    /* Here are the known values.  Note that these get written to
     * /gammaUpdowndown because I can't figure out how to make Petsc create a
     * new dataset, or a group, or anything I can set attributes on. */
    
    WRITE_PARAM_INT(COMPUTE_DIM);
    WRITE_PARAM_INT(N1);
    WRITE_PARAM_INT(N2);
    WRITE_PARAM_DOUBLE(DT);
    WRITE_PARAM_DOUBLE(START_TIME);
    WRITE_PARAM_DOUBLE(FINAL_TIME);
    WRITE_PARAM_INT(RHO);
    WRITE_PARAM_INT(UU);
    WRITE_PARAM_INT(U1);
    WRITE_PARAM_INT(U2);
    WRITE_PARAM_INT(U3);
    WRITE_PARAM_INT(B1);
    WRITE_PARAM_INT(B2);
    WRITE_PARAM_INT(B3);

    writeProblemSpecificData(parametersViewer, ts->problemSpecificData);

    PetscViewerDestroy(&parametersViewer);

    VecDestroy(&gCovPetscVec);
    VecDestroy(&gConPetscVec);

    DMDestroy(&metricDMDA);

    PetscPrintf(PETSC_COMM_WORLD, "done\n");
  }

  if (ts->t > ts->tDump || fabs(ts->t - ts->tDump) < 1e-10
                        || fabs(ts->t - FINAL_TIME) < 1e-10 
     )
  {
    char primVarsFileName[50], residualsFileName[50];
    sprintf(primVarsFileName, "%s%06d.h5", DUMP_FILE_PREFIX, ts->dumpCounter);
    sprintf(residualsFileName, "%s%06d.h5", RESIDUALS_DUMP_FILE_PREFIX,
            ts->dumpCounter
           );

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, primVarsFileName,
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) ts->primPetscVec, "primVars");
    VecView(ts->primPetscVec, viewer);
    PetscViewerDestroy(&viewer);

    /* Get the residual at the last iteration of the SNES solver */
    Vec residualPetscVec;
    SNESGetFunction(ts->snes, &residualPetscVec, NULL, NULL);

    PetscViewerHDF5Open(PETSC_COMM_WORLD, residualsFileName,
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) residualPetscVec, "residuals");
    VecView(residualPetscVec, viewer);
    PetscViewerDestroy(&viewer);

    PetscPrintf(PETSC_COMM_WORLD,
                "\nDumped %s and %s at t = %f\n\n",
                ts->t, primVarsFileName, residualsFileName);
    
    ts->tDump += DT_DUMP;
    ts->dumpCounter++;
  }
  
}
