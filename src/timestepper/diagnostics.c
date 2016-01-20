#include "timestepper.h"

void diagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  if (ts->timeStepCounter==0)
  {
    /* Dump geometry data */

    PetscPrintf(PETSC_COMM_WORLD, 
                "Dumping the grid data..."
               );

    DM metricDM;
    Vec gCovVec, gConVec;

    DM gDetDM, alphaDM;
    Vec gDetVec, alphaVec;

    DM coordinatesDM;
    Vec xCoordsVec, XCoordsVec, cartesianCoordsVec;

#if (COMPUTE_DIM==1)

    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, 16, 0, NULL,
                 &metricDM
                );
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, 1, 0, NULL,
                 &gDetDM
                );
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, 1, 0, NULL,
                 &alphaDM
                );
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N1, 4, 0, NULL,
                 &coordinatesDM
                );

#elif (COMPUTE_DIM==2)

    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX, N1, N2, PETSC_DECIDE, PETSC_DECIDE,
                 16, 0, PETSC_NULL, PETSC_NULL, &metricDM
                );
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX, N1, N2, PETSC_DECIDE, PETSC_DECIDE,
                 1, 0, PETSC_NULL, PETSC_NULL, &gDetDM
                );
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX, N1, N2, PETSC_DECIDE, PETSC_DECIDE,
                 1, 0, PETSC_NULL, PETSC_NULL, &alphaDM
                );
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX, N1, N2, PETSC_DECIDE, PETSC_DECIDE,
                 4, 0, PETSC_NULL, PETSC_NULL, &coordinatesDM
                );

#endif /* Choose dimension */

    DMCreateGlobalVector(metricDM, &gCovVec);
    DMCreateGlobalVector(metricDM, &gConVec);

    DMCreateGlobalVector(gDetDM,  &gDetVec);
    DMCreateGlobalVector(alphaDM, &alphaVec);

    DMCreateGlobalVector(coordinatesDM, &xCoordsVec);
    DMCreateGlobalVector(coordinatesDM, &XCoordsVec);
    DMCreateGlobalVector(coordinatesDM, &cartesianCoordsVec);

    ARRAY(gCovGlobal);
    ARRAY(gConGlobal);
    DMDAVecGetArrayDOF(metricDM, gCovVec, &gCovGlobal);
    DMDAVecGetArrayDOF(metricDM, gConVec, &gConGlobal);

    ARRAY(gDetGlobal);
    ARRAY(alphaGlobal);
    DMDAVecGetArrayDOF(gDetDM,  gDetVec,  &gDetGlobal);
    DMDAVecGetArrayDOF(alphaDM, alphaVec, &alphaGlobal);

    ARRAY(xCoordsGlobal);
    ARRAY(XCoordsGlobal);
    ARRAY(cartesianCoordsGlobal);
    DMDAVecGetArrayDOF(coordinatesDM, xCoordsVec, &xCoordsGlobal);
    DMDAVecGetArrayDOF(coordinatesDM, XCoordsVec, &XCoordsGlobal);
    DMDAVecGetArrayDOF(coordinatesDM, cartesianCoordsVec, 
                       &cartesianCoordsGlobal
                      );

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

        REAL XCoords[NDIM], xCoords[NDIM], cartesianCoords[NDIM];
        getXCoords(&zone, CENTER, XCoords);
        XTox(XCoords, xCoords);
        XToCartesian(XCoords, cartesianCoords);

        struct geometry geom;
        setGeometry(XCoords, &geom);

        INDEX_PETSC(gDetGlobal, &zone, 0) = geom.gDet;

        for (int mu=0; mu<NDIM; mu++)
        {
          INDEX_PETSC(xCoordsGlobal, &zone, mu) = xCoords[mu];
          INDEX_PETSC(XCoordsGlobal, &zone, mu) = XCoords[mu];
          INDEX_PETSC(cartesianCoordsGlobal, &zone, mu) =
            cartesianCoords[mu];
          
          for (int nu=0; nu<NDIM; nu++)
          {
            INDEX_PETSC(gCovGlobal, &zone, nu + NDIM*mu) = geom.gCov[mu][nu];
            INDEX_PETSC(gConGlobal, &zone, nu + NDIM*mu) = geom.gCon[mu][nu];
          }
        }

      }
    }

    DMDAVecRestoreArrayDOF(coordinatesDM, xCoordsVec, &xCoordsGlobal);
    DMDAVecRestoreArrayDOF(coordinatesDM, XCoordsVec, &XCoordsGlobal);
    DMDAVecRestoreArrayDOF(coordinatesDM, cartesianCoordsVec, 
                           &cartesianCoordsGlobal
                          );

    DMDAVecRestoreArrayDOF(metricDM, gCovVec, &gCovGlobal);
    DMDAVecRestoreArrayDOF(metricDM, gConVec, &gConGlobal);
    DMDAVecRestoreArrayDOF(gDetDM,  gDetVec,  &gDetGlobal);
    DMDAVecRestoreArrayDOF(alphaDM, alphaVec, &alphaGlobal);

    char gridFilename[50];
    sprintf(gridFilename, "grid.h5");

    PetscViewer viewer;
    
    PetscViewerHDF5Open(PETSC_COMM_WORLD, gridFilename,
                        FILE_MODE_WRITE, &viewer);

    PetscObjectSetName((PetscObject) gCovVec, "gCov");
    PetscObjectSetName((PetscObject) gConVec, "gCon");
    PetscObjectSetName((PetscObject) ts->connectionPetscVec,
                       "gammaUpdowndown");

    PetscObjectSetName((PetscObject) gDetVec,  "gDet");
    PetscObjectSetName((PetscObject) alphaVec, "alphaLapse");

    PetscObjectSetName((PetscObject) xCoordsVec, "xCoords");
    PetscObjectSetName((PetscObject) XCoordsVec, "XCoords");
    PetscObjectSetName((PetscObject) cartesianCoordsVec, "cartesianCoords");

    VecView(gCovVec,                viewer);
    VecView(gConVec,                viewer);
    VecView(ts->connectionPetscVec, viewer);
    VecView(gDetVec,                viewer);
    VecView(alphaVec,               viewer);
    VecView(xCoordsVec,             viewer);
    VecView(XCoordsVec,             viewer);
    VecView(cartesianCoordsVec,     viewer);

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

    PetscViewerHDF5PushGroup(viewer, "/cartesianCoords");
    WRITE_PARAM_INT(COMPUTE_DIM,    viewer);
    WRITE_PARAM_INT(N1,             viewer);
    WRITE_PARAM_INT(N2,             viewer);
    WRITE_PARAM_DOUBLE(DT,          viewer);
    WRITE_PARAM_DOUBLE(START_TIME,  viewer);
    WRITE_PARAM_DOUBLE(FINAL_TIME,  viewer);
    WRITE_PARAM_INT(RHO,            viewer);
    WRITE_PARAM_INT(UU,             viewer);
    WRITE_PARAM_INT(U1,             viewer);
    WRITE_PARAM_INT(U2,             viewer);
    WRITE_PARAM_INT(U3,             viewer);
    WRITE_PARAM_INT(B1,             viewer);
    WRITE_PARAM_INT(B2,             viewer);
    WRITE_PARAM_INT(B3,             viewer);
    #if (CONDUCTION)
      WRITE_PARAM_INT(PHI,          viewer);
    #endif
    #if (VISCOSITY)
      WRITE_PARAM_INT(PSI,          viewer);
    #endif
    WRITE_PARAM_INT(DOF,            viewer);

    writeProblemSpecificData(viewer, ts->problemSpecificData);
    PetscViewerHDF5PopGroup(viewer);

    PetscViewerDestroy(&viewer);

    VecDestroy(&gCovVec);
    VecDestroy(&gConVec);
    VecDestroy(&gDetVec);
    VecDestroy(&alphaVec);

    VecDestroy(&xCoordsVec);
    VecDestroy(&XCoordsVec);
    VecDestroy(&cartesianCoordsVec);

    DMDestroy(&metricDM);
    DMDestroy(&gDetDM);
    DMDestroy(&alphaDM);
    DMDestroy(&coordinatesDM);

    PetscPrintf(PETSC_COMM_WORLD, "done\n");
  }

  if (ts->t > ts->tDump || fabs(ts->t - ts->tDump) < 1e-10
                        || fabs(ts->t - FINAL_TIME) < 1e-10 
     )
  {
    char dumpFileName[50];
    sprintf(dumpFileName, "%s%06d.h5", DUMP_FILE_PREFIX, ts->dumpCounter);

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, dumpFileName,
                        FILE_MODE_WRITE, &viewer
                       );

    /* Output the solved primitive variables */
    PetscObjectSetName((PetscObject) ts->primPetscVec, "primVars");
    VecView(ts->primPetscVec, viewer);

    /* Output the old primitive variables. Need this to compute time
     * derivatives */
    PetscObjectSetName((PetscObject) ts->primPetscVecOld, "primVarsOld");
    VecView(ts->primPetscVecOld, viewer);

    /* Get the residual at the last iteration of the SNES solver and output that */
    Vec residualPetscVec;
    SNESGetFunction(ts->snes, &residualPetscVec, NULL, NULL);
    PetscObjectSetName((PetscObject) residualPetscVec, "residuals");
    VecView(residualPetscVec, viewer);

    PetscViewerHDF5PushGroup(viewer, "/primVars");
    REAL dt = ts->dt; WRITE_PARAM_DOUBLE(dt, viewer);
    REAL t  = ts->t;  WRITE_PARAM_DOUBLE(t, viewer);

    writeProblemSpecificData(viewer, ts->problemSpecificData);
    PetscViewerHDF5PopGroup(viewer);

    PetscViewerDestroy(&viewer);

    PetscPrintf(PETSC_COMM_WORLD,
                "\nDumped %s at t = %f\n\n",
                dumpFileName, ts->t
               );
    
    ts->tDump += DT_DUMP;
    ts->dumpCounter++;
  }
  
}
