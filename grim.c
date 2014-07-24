#include "grim.h"

int main(int argc, char **argv)
{
    TS ts;
    SNES snes;
    Vec soln;
    DM dmda;

    int X1Start, X2Start;
    int X1Size, X2Size;

    PetscInitialize(&argc, &argv, PETSC_NULL, help);

    DMDACreate2d(PETSC_COMM_WORLD, 
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DMDA_STENCIL_STAR,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 DOF, NG, PETSC_NULL, PETSC_NULL, &dmda);

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    SetCoordinates(dmda);
//    DMDASetUniformCoordinates(dmda, 
//                              DX1/2., 1. - DX1/2.,
//                              DX2/2., 1. - DX2/2., 
//                              0., 0.);


    DMDASetFieldName(dmda, RHO, "Density");
    DMDASetFieldName(dmda, UU, "Internal Energy");
    DMDASetFieldName(dmda, U1, "Ur");
    DMDASetFieldName(dmda, U2, "Utheta");
    DMDASetFieldName(dmda, U3, "Uphi");
    DMDASetFieldName(dmda, B1, "Br");
    DMDASetFieldName(dmda, B2, "Btheta");
    DMDASetFieldName(dmda, B3, "Bphi");    

    DMCreateGlobalVector(dmda, &soln);

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetDM(ts, dmda);
    TSSetIFunction(ts, PETSC_NULL, ComputeResidual, NULL);


    InitialCondition(ts, soln);

    TSSetSolution(ts, soln);
    TSMonitorSet(ts, Monitor, NULL, NULL);
    TSGetSNES(ts, &snes);
    SNESMonitorSet(snes, SNESMonitor, NULL, NULL);
    TSSetType(ts, TSTHETA);
    TSSetFromOptions(ts);

    TSSolve(ts, soln);

    TSDestroy(&ts);
    VecDestroy(&soln);
    DMDestroy(&dmda);

    PetscFinalize();
    return(0);
}

PetscErrorCode ComputeResidual(TS ts,
                               PetscScalar t,
                               Vec primPetscVec, Vec dPrimDtPetscVec,
                               Vec residualPetscVec, void *ptr)
{
    DM da;
    Vec localPrimPetscVec;
    int X1Start, X2Start, X3Start;
    int X1Size, X2Size, X3Size;
    int X1End, X2End, X3End;

    /* localPrimPetscVec is padded with ghost zones */
    DMGetLocalVector(da, &localPrimPetscVec);
    
    DMDAGetCorners(da, &X1Start, &X2Start, &X3Start,
                       &X1Size, &X2Size, &X3Size);

    X1End = X1Start + X1Size;
    X2End = X2Start + X2Size;
    X3End = X3Start + X3Size;

    /* Only transfer ghost zone information between MPI nodes for primPetscVec
     * into localPrimPetscVec. dPrimDtPetscVec and residualPetscVec are global
     * vecs and do not have ghost zones. Ghost zone info not needed for
     * dPrimDtPetscVec and residualPetscVec */
    DMGlobalToLocalBegin(da, primPetscVec, INSERT_VALUES, localPrimPetscVec);
    DMGlobalToLocalEnd(da, primPetscVec, INSERT_VALUES, localPrimPetscVec);

    /* Get the data on the local node as a pointer. In the following,  localPrim
     * is padded with ghost zones whereas localDPrimDt and localResidual are
     * not. localDPrimDt and localResidual get the data of the global vectors
     * dPrimDtPetscVec and residualPetscVec which lie on the current node. */
    REAL *localPrim, *localDPrimDt, *localResidual;
    VecGetArray(localPrimPetscVec, &localPrim);
    VecGetArray(dPrimDtPetscVec, &localDPrimDt);
    VecGetArray(residualPetscVec, &localResidual);

    int numTilesInX2 = X2End/TILE_SIZE_X2;
    int numTilesInX1 = X1End/TILE_SIZE_X1;

    for (int tileNumX2=0; tileNumX2<numTilesInX2; tileNumX2++)
    {
        for (int tileNumX1=0; tileNumX1<numTilesInX1; tileNumX1++)
        {
            struct gridZone2D zone;

            REAL primTile[ (TILE_SIZE_X1 + 2*NG)\
                          *(TILE_SIZE_X2 + 2*NG)*DOF];

            for (zone.jTile=0; zone.jTile<TILE_SIZE_X2; zone.jTile++)
            {
                for (zone.iTile=0; zone.iTile<TILE_SIZE_X1; zone.iTile++)
                {
                    zone.i = X1Start + iTile + tileNumX1*TILE_SIZE_X1;
                    zone.j = X2Start + jTile + tileNumX2*TILE_SIZE_X2;

                    for (int var=0; var<DOF; var++)
                    {
                        primTile[INDEX_LOCAL(zone.iTile, zone.jTile, var)] =
                        localPrim[INDEX_GLOBAL(zone.i, zone.j, var)];
                    }

                }
            }

            for (zone.jTile=0; zone.jTile<TILE_SIZE_X2; zone.jTile++)
            {
                for (zone.iTile=0; zone.iTile<TILE_SIZE_X1; zone.iTile++)
                {
                    zone.i = X1Start + iTile + tileNumX1*TILE_SIZE_X1;
                    zone.j = X2Start + jTile + tileNumX2*TILE_SIZE_X2;

                    setZone2DBoundaryFlags(zone);

                    applyTileBoundaryConditions(zone, localPrim, primTile);
                }
            }

            for (zone.jTile=0; zone.jTile<TILE_SIZE_X2; zone.jTile++)
            {
                for (zone.iTile=0; zone.iTile<TILE_SIZE_X1; zone.iTile++)
                {


                }
            }

        }
    }




    VecRestoreArray(primPetscVec, &prim);
    VecRestoreArray(dPrimDtPetscVec, &dprimDt);
    VecRestoreArray(residualPetscVec, &residual);

    DMRestoreLocalVector(da, &localPrimPetscVec);

    return(0);








    PetscScalar *prim, *dprim_dt, *f;
    VecGetArray(Prim, &prim);
    VecGetArray(dPrim_dt, &dprim_dt);
    VecGetArray(F, &f);

    cl::Buffer primBuffer, dprimBuffer_dt, fbuffer;
    PetscInt size = DOF*N1*N2*sizeof(PetscScalar);

    primBuffer = cl::Buffer(context,
                            CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                            size, &(prim[0]), &clErr);
    dprimBuffer_dt = cl::Buffer(context,
                                CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                size, &(dprim_dt[0]), &clErr);
    fbuffer = cl::Buffer(context,
                         CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
                         size, &(f[0]), &clErr);


    clErr = kernel.setArg(0, primBuffer);
    clErr = kernel.setArg(1, dprimBuffer_dt);
    clErr = kernel.setArg(2, fbuffer);

    cl::NDRange global(N1, N2);
    cl::NDRange local(TILE_SIZE_X1, TILE_SIZE_X2);
    clErr = queue.enqueueNDRangeKernel(kernel,
                                       cl::NullRange,
                                       global, local,
                                       NULL, NULL);

    f = (PetscScalar*)queue.enqueueMapBuffer(fbuffer,
                                             CL_FALSE,
                                             CL_MAP_READ,
                                             0, size,
                                             NULL, NULL, &clErr);

    clErr = queue.finish();

    VecRestoreArray(Prim, &prim);
    VecRestoreArray(dPrim_dt, &dprim_dt);
    VecRestoreArray(F, &f);

    return(0.);
}


//PetscErrorCode Monitor(TS ts, 
//                       PetscInt step,
//                       PetscReal t,
//                       Vec Prim,
//                       void *ptr)
//{
//    DM dmda;
//    TSGetDM(ts, &dmda);
//    PetscScalar ***prim;
//    PetscScalar gcov[NDIM][NDIM];
//
//    DMDAVecGetArrayDOF(dmda, Prim, &prim);
//
//    for (int j=0; j<N2; j++) {
//        for (int i=0; i<N1; i++) {        
//
//            REAL X1 = i_TO_X1_CENTER(i);
//            REAL X2 = j_TO_X2_CENTER(j);
//            gCovCalc(gcov, X1, X2);
//
//            REAL vars[DOF];
//            for (int var=0; var<DOF; var++)
//                vars[var] = prim[j][i][var];
//
//            REAL gamma;
//            gammaCalc(&gamma, vars, gcov);
//
//            setFloorInit(vars, gamma, X1, X2);
//
//            for (int var=0; var<DOF; var++)
//                prim[j][i][var] = vars[var];
//        }
//    }
//
//    DMDAVecRestoreArrayDOF(dmda, Prim, &prim);
//
//    REAL dt, dtDump;
//    TSGetTimeStep(ts, &dt);
//    
//    static PetscInt counter = 0;
//    static PetscScalar tDump = 0.;
//    dtDump = 1.;
//
//    SNES snes;
//    TSGetSNES(ts, &snes);
//    Vec Res;
//    SNESGetFunction(snes, &Res, NULL, NULL);
//
//    if (t > tDump) {
//        printf("Dumping data..\n");
//        char filename[50];
//        char errname[50];
//        sprintf(filename, "plot%04d.h5", counter);
//        sprintf(errname, "residual%04d.h5", counter);
//
//        PetscViewer viewer;
//        PetscViewerHDF5Open(PETSC_COMM_WORLD, filename,
//                            FILE_MODE_WRITE, &viewer);
//        PetscObjectSetName((PetscObject) Prim, "soln");
//        VecView(Prim, viewer);
//        PetscViewerDestroy(&viewer);
//
//        PetscViewerHDF5Open(PETSC_COMM_WORLD, errname,
//                            FILE_MODE_WRITE, &viewer);
//        PetscObjectSetName((PetscObject) Res, "Err");
//        VecView(Res, viewer);
//        PetscViewerDestroy(&viewer);
//
//        tDump = tDump + dtDump;
//        counter++;
//    }
//
//
//    TSSetTimeStep(ts, .03);
//
//    return(0.);
//}
//
//PetscErrorCode SNESMonitor(SNES snes, PetscInt its, PetscReal norm, void *ptr)
//{
//    if (its==0) {
//        Vec Res;
//        SNESGetFunction(snes, &Res, NULL, NULL);
//
//        char errname[50];
//        sprintf(errname, "SNESresidual.h5");
//
//        PetscViewer viewer;
//        PetscViewerHDF5Open(PETSC_COMM_WORLD, errname,
//                            FILE_MODE_WRITE, &viewer);
//        PetscObjectSetName((PetscObject) Res, "Err");
//        VecView(Res, viewer);
//        PetscViewerDestroy(&viewer);
//    }
//
//    return(0);
//}
