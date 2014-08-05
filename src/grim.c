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
