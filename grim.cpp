#include "grim.h"

PetscErrorCode SetCoordinates(DM dmda)
{
    DM coordDM;
    Vec coordVec;

    int X1Start, X2Start;
    int X1Size, X2Size;

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    DMDASetUniformCoordinates(dmda,
                              0., 1., 0., 1., 0., 1.);
    DMGetCoordinateDM(dmda, &coordDM);
    DMGetCoordinatesLocal(dmda, &coordVec);

    PetscScalar X1, X2, r, theta;
    DMDACoor2d **coord;
    
    DMDAVecGetArray(coordDM, coordVec, &coord);

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
    
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            BLCoords(&r, &theta, X1, X2);

            coord[j][i].x = r*PetscSinScalar(theta);
            coord[j][i].y = r*PetscCosScalar(theta);
        }
    }
    DMDAVecRestoreArray(coordDM, coordVec, &coord);

    DMSetCoordinatesLocal(dmda, coordVec);

    return (0);
}

int main(int argc, char **argv)
{
    TS ts;
    Vec soln;
    DM dmda;

    int X1Start, X2Start;
    int X1Size, X2Size;

    PetscInitialize(&argc, &argv, PETSC_NULL, help);

//    viennacl::ocl::set_context_platform_index(0, 1);

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

//    std::ifstream sourceFile("computeresidual.cl");
//    std::string sourceCode((std::istreambuf_iterator<char>(sourceFile)),
//                            std::istreambuf_iterator<char>());
//
//    std::string BuildOptions("\
//                              -D X1_SIZE=" +
//                             std::to_string(X1Size) +
//                             " -D X2_SIZE=" + 
//                             std::to_string(X2Size) +
//                             " -D TOTAL_X1_SIZE=" + 
//                             std::to_string(X1Size+2*NG) + 
//                             " -D TOTAL_X2_SIZE=" +
//                             std::to_string(X2Size+2*NG));
//
//    viennacl::ocl::current_context().build_options(BuildOptions);
//
//    PetscScalar start = std::clock();
//    program =
//        viennacl::ocl::current_context().add_program(sourceCode,
//                                                     "computeresidual");
//    PetscScalar end = std::clock();
//    PetscScalar time = (end - start)/(PetscScalar)CLOCKS_PER_SEC;
//    PetscPrintf(PETSC_COMM_WORLD, 
//                "Time taken for kernel compilation = %f\n", time);
//    
//    program.add_kernel("ComputeResidual");

    clErr = cl::Platform::get(&platforms);
    CheckCLErrors(clErr, "cl::Platform::get");

    clErr = platforms.at(1).getDevices(CL_DEVICE_TYPE_CPU, &devices);
    CheckCLErrors(clErr, "cl::Platform::getDevices");

    context = cl::Context(devices, NULL, NULL, NULL, &clErr);
    CheckCLErrors(clErr, "cl::Context::Context");

    queue = cl::CommandQueue(context, devices.at(0), 0, &clErr);
    CheckCLErrors(clErr, "cl::CommandQueue::CommandQueue");

    std::ifstream sourceFile("computeresidual.cl");
    std::string sourceCode((std::istreambuf_iterator<char>(sourceFile)),
                            std::istreambuf_iterator<char>());
    cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(),
                                sourceCode.length()+1));
    
    program = cl::Program(context, source, &clErr);
    CheckCLErrors(clErr, "cl::Program::Program");

    std::string BuildOptions("\
                              -D X1_SIZE=" +
                             std::to_string(X1Size) +
                             " -D X2_SIZE=" + 
                             std::to_string(X2Size) +
                             " -D TOTAL_X1_SIZE=" + 
                             std::to_string(X1Size+2*NG) + 
                             " -D TOTAL_X2_SIZE=" +
                             std::to_string(X2Size+2*NG));

    PetscScalar start = std::clock();
    clErr = program.build(devices, BuildOptions.c_str(), NULL, NULL);
    const char *buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(
                                                devices.at(0),
                                                &clErr).c_str();
    PetscPrintf(PETSC_COMM_WORLD, "%s\n", buildlog);
    CheckCLErrors(clErr, "cl::Program::build");
    PetscScalar end = std::clock();

    PetscScalar time = (end - start)/(PetscScalar)CLOCKS_PER_SEC;
    PetscPrintf(PETSC_COMM_WORLD, 
                "Time taken for kernel compilation = %f\n", time);

    kernel = cl::Kernel(program, "ComputeResidual", &clErr);
    CheckCLErrors(clErr, "cl::Kernel::Kernel");

    cl_ulong localMemSize = kernel.getWorkGroupInfo<CL_KERNEL_LOCAL_MEM_SIZE>(
                                        devices.at(0), &clErr);
    cl_ulong privateMemSize = kernel.getWorkGroupInfo<CL_KERNEL_PRIVATE_MEM_SIZE>(
                                        devices.at(0), &clErr);
    printf("Local memory used = %llu\n", (unsigned long long)localMemSize);
    printf("Private memory used = %llu\n", (unsigned long long)privateMemSize);

//    PetscViewer viewer;
//    PetscViewerHDF5Open(PETSC_COMM_WORLD,"init.h5",
//                        FILE_MODE_READ, &viewer);
//    PetscObjectSetName((PetscObject) soln,"soln");
//    VecLoad(soln, viewer);
//    PetscViewerDestroy(&viewer);

    InitialCondition(ts, soln);

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD,"plot0.h5",
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) soln, "soln");
    VecView(soln, viewer);
    PetscViewerDestroy(&viewer);


//    Vec X;
//    VecDuplicate(soln, &X);
//    VecSet(X, 10.);
//    PetscScalar ***x, ***prim;
//    DMDAVecGetArrayDOF(dmda, X, &x);
//    DMDAVecGetArrayDOF(dmda, soln, &prim);
//    for (int j=0; j<N2; j++)
//        for (int i=0; i<N1; i++) {
//            PetscScalar X1 = i_TO_X1_CENTER(i);
//            PetscScalar X2 = j_TO_X2_CENTER(j);
//
//            PetscScalar sources[DOF], var[DOF];
//            for (int n=0; n<DOF; n++) {
//                sources[n] = 0.;
//                var[n] = prim[j][i][n];
//            }
//            
//            PetscScalar gcon[NDIM][NDIM], gcov[NDIM][NDIM], gdet;
//            PetscScalar ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
//            PetscScalar mhd[NDIM][NDIM], flux[DOF], U[DOF];
//            PetscScalar alpha, gamma, g;
//            PetscScalar cmin, cmax, bsqr;
//
//            gCovCalc(gcov, X1, X2);
//            gDetCalc(&gdet, gcov);
//            gConCalc(gcon, gcov, gdet);
//            g = sqrt(-gdet);
//
//            gammaCalc(&gamma, var, gcov);
//            alphaCalc(&alpha, gcon);
//            uconCalc(ucon, gamma, alpha, var, gcon);
//            covFromCon(ucov, ucon, gcov);
//            bconCalc(bcon, var, ucon, ucov);
//            covFromCon(bcov, bcon, gcov);
//
//            mhdCalc(mhd, var, ucon, ucov, bcon, bcov);
//
//            addSources(sources, 
//                       ucon, ucov, bcon, bcov, 
//                       gcon, gcov, mhd, var, g,
//                       X1, X2);
//
//
//            ComputeFluxAndU(flux, U,
//                            ucon, ucov,
//                            bcon, bcov,
//                            gcon, gcov,
//                            mhd, var, g, 1);
//
//            conDotCov(&bsqr, bcon, bcov);
//            VChar(&cmin, &cmax, ucon, ucov, bsqr, gcon, var, 1); 
//
//            x[j][i][RHO] = cmin;
//            x[j][i][UU] = cmax;
//            x[j][i][U1] = flux[U1];
//            x[j][i][U2] = flux[U2];
//            x[j][i][U3] = flux[U3];
//            x[j][i][B1] = flux[B1];
//            x[j][i][B2] = flux[B2];
//            x[j][i][B3] = flux[B3];
//
//        }
//    DMDAVecRestoreArrayDOF(dmda, X, &x);
//    DMDAVecRestoreArrayDOF(dmda, soln, &prim);
//
//    PetscViewer viewer;
//    PetscViewerHDF5Open(PETSC_COMM_WORLD,"metric.h5",
//                        FILE_MODE_WRITE, &viewer);
//    PetscObjectSetName((PetscObject) X, "gcov00");
//    VecView(X, viewer);
//    PetscViewerDestroy(&viewer);
//    VecDestroy(&X);
        

//    Benchmark(ts, soln);

    TSSetSolution(ts, soln);
    TSMonitorSet(ts, Monitor, NULL, NULL);
    TSSetType(ts, TSTHETA);
    TSSetFromOptions(ts);

    TSSolve(ts, soln);

    TSDestroy(&ts);
    VecDestroy(&soln);
    DMDestroy(&dmda);

    PetscFinalize();
    return(0);
}

//PetscErrorCode ComputeResidual(TS ts,
//                               PetscScalar t,
//                               Vec Prim, Vec dPrim_dt,
//                               Vec F, void *ptr)
//{
//    const viennacl::vector<PetscScalar> *prim, *dprim_dt;
//    viennacl::vector<PetscScalar> *f;
//
//    VecViennaCLGetArrayRead(Prim, &prim);
//    VecViennaCLGetArrayRead(dPrim_dt, &dprim_dt);
//    VecViennaCLGetArrayWrite(F, &f);
//
//    viennacl::ocl::kernel kernel = program.get_kernel("ComputeResidual");
//
//    kernel.local_work_size(0, TILE_SIZE_X1);
//    kernel.local_work_size(1, TILE_SIZE_X2);
//    kernel.global_work_size(0, N1);
//    kernel.global_work_size(1, N2);
//
//    viennacl::ocl::enqueue(kernel(*prim, *dprim_dt, *f));
//    viennacl::backend::finish();
//
//    VecViennaCLRestoreArrayRead(Prim, &prim);
//    VecViennaCLRestoreArrayRead(dPrim_dt, &dprim_dt);
//    VecViennaCLRestoreArrayWrite(F, &f);
//
//    return (0);
//}

PetscErrorCode ComputeResidual(TS ts,
                               PetscScalar t,
                               Vec Prim, Vec dPrim_dt,
                               Vec F, void *ptr)
{
//    DM dmda;
//    int X1Start, X1Size;
//    int X2Start, X2Size;
//    TSGetDM(ts, &dmda);
//
//    DMDAGetCorners(dmda, 
//                   &X1Start, &X2Start, NULL,
//                   &X1Size, &X2Size, NULL);
//
//    Vec PrimLocal;
//    DMGetLocalVector(dmda, &PrimLocal);
//
//    DMGlobalToLocalBegin(dmda, Prim, INSERT_VALUES, PrimLocal);
//    DMGlobalToLocalEnd(dmda, Prim, INSERT_VALUES, PrimLocal);
//
//
//    PetscScalar ***prim, ***dprim_dt, ***f;
//    DMDAVecGetArrayDOF(dmda, PrimLocal, &prim);
//    DMDAVecGetArrayDOF(dmda, dPrim_dt, &dprim_dt);
//    DMDAVecGetArrayDOF(dmda, F, &f);
//
//    cl::Buffer primBuffer, dprimBuffer_dt, fbuffer;
//    PetscInt sizeWithoutNG = DOF*N1*N2*sizeof(PetscScalar);
//    PetscInt sizeWithNG = DOF*(N1+2*NG)*(N2+2*NG)*sizeof(PetscScalar);
//
//    primBuffer = cl::Buffer(context,
//                            CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
//                            sizeWithNG, &(prim[-2][-2][0]), &clErr);
//    dprimBuffer_dt = cl::Buffer(context,
//                                CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
//                                sizeWithoutNG, &(dprim_dt[0][0][0]), &clErr);
//    fbuffer = cl::Buffer(context,
//                         CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
//                         sizeWithoutNG, &(f[0][0][0]), &clErr);
//
//
//    clErr = kernel.setArg(0, primBuffer);
//    clErr = kernel.setArg(1, dprimBuffer_dt);
//    clErr = kernel.setArg(2, fbuffer);
//
//    cl::NDRange global(N1, N2);
//    cl::NDRange local(TILE_SIZE_X1, TILE_SIZE_X2);
//    clErr = queue.enqueueNDRangeKernel(kernel,
//                                       cl::NullRange,
//                                       global, local,
//                                       NULL, NULL);
//
////    f = (PetscScalar***)queue.enqueueMapBuffer(fbuffer,
////                                             CL_FALSE,
////                                             CL_MAP_READ,
////                                             0, size,
////                                             NULL, NULL, &clErr);
//
//    clErr = queue.finish();
//
//    DMDAVecRestoreArrayDOF(dmda, PrimLocal, &prim);
//    DMDAVecRestoreArrayDOF(dmda, dPrim_dt, &dprim_dt);
//    DMDAVecRestoreArrayDOF(dmda, F, &f);
//
//    DMRestoreLocalVector(dmda, &PrimLocal);


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

//    DM dmda;
//    int X1Start, X1Size;
//    int X2Start, X2Size;
//    TSGetDM(ts, &dmda);
//
//    DMDAGetCorners(dmda, 
//                   &X1Start, &X2Start, NULL,
//                   &X1Size, &X2Size, NULL);
//
//    Vec PrimLocal;
//    DMGetLocalVector(dmda, &PrimLocal);
//
//    DMGlobalToLocalBegin(dmda, Prim, INSERT_VALUES, PrimLocal);
//    DMGlobalToLocalEnd(dmda, Prim, INSERT_VALUES, PrimLocal);
//
//    PetscScalar ***prim, ***dprim_dt, ***f;
//    DMDAVecGetArrayDOF(dmda, PrimLocal, &prim);
//    DMDAVecGetArrayDOF(dmda, dPrim_dt, &dprim_dt);
//    DMDAVecGetArrayDOF(dmda, F, &f);
//
//    REAL primTile[(TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*DOF];
//
//    REAL dU_dt[DOF], primEdge[DOF];
//    REAL fluxL[DOF], fluxR[DOF];
//    REAL uL[DOF], uR[DOF];
//    REAL fluxX1L[DOF], fluxX1R[DOF];
//    REAL fluxX2L[DOF], fluxX2R[DOF];
//    REAL X1, X2;
//
//    // Geometry
//    REAL gcon[NDIM][NDIM], gcov[NDIM][NDIM];
//
//    //Physics
//    REAL ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
//    REAL mhd[NDIM][NDIM];
//    REAL vars[DOF], dvars_dt[DOF];
//
//    for (int j=0; j<N2; j++)
//        for (int i=0; i<N1; i++) {
//            
//            for (int var=0; var<DOF; var++) {
//                vars[var] = prim[j][i][DOF];
//                dU_dt[var] = 0.;
//            }
//
//            X1 = i_TO_X1_FACE(i); X2 = j_TO_X2_CENTER(j);
//            addSources(dU_dt,
//                       ucon, ucov, bcon, bcov,
//                       gcon, gcov, mhd, vars, 
//                       X1, X2);
//
//            for (int var=0; var<DOF; var++)
//                f[j][i][var] = dU_dt[var];
//        }
//
////    for (int j=0; j<N2; j++) {
////        for (int iNg=-NG; iNg<0; iNg++) {
////            for (int var=0; var<DOF; var++) {
////                prim[j][iNg][var] = prim[j][0][var];
////            }
////        }
////
////        for (int iNg=N1; iNg<N1+NG; iNg++) {
////            for (int var=0; var<DOF; var++) {
////                prim[j][iNg][var] = prim[j][N1-1][var];
////            }
////        }
////    }
////
////    for (int jNg=-NG; jNg<0; jNg++) {
////        for (int i=0; i<N1; i++) {
////            for (int var=0; var<DOF; var++) {
////                prim[jNg][i][var] = prim[-jNg-1][i][var];
////            }
////        }
////    }
////
////    for (int jNg=N1; jNg<N2+NG; jNg++) {
////        for (int i=0; i<N1; i++) {
////            for (int var=0; var<DOF; var++) {
////                prim[jNg][i][var] = prim[-jNg+2*N1-1][i][var];
////            }
////        }
////    }
////
////    for (int jB=0; jB<N2/TILE_SIZE_X2; jB++) {
////        for (int iB=0; iB<N1/TILE_SIZE_X1; iB++) {
////            for (int jTile=-NG; jTile<TILE_SIZE_X2+NG; jTile++) {
////                for (int iTile=-NG; iTile<TILE_SIZE_X1+NG; iTile++) {
////
////                        int i = iTile + iB*TILE_SIZE_X1;
////                        int j = jTile + jB*TILE_SIZE_X2;
////
////                    for (int var=0; var<DOF; var++) {
////                        primTile[INDEX_LOCAL(iTile,jTile,var)] = 
////                        prim[j][i][var];
////                    }
////                }
////            }
////
////            for (int jTile=0; jTile<TILE_SIZE_X2; jTile++) {
////                for (int iTile=0; iTile<TILE_SIZE_X1; iTile++) {
////
////                    int i = iTile + iB*TILE_SIZE_X1;
////                    int j = jTile + jB*TILE_SIZE_X2;
////
////                    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j);
////                    //ComputedU_dt
////
////                    for (int var=0; var<DOF; var++) {
////                        dU_dt[var] = 0.;
////                        vars[var] = primTile[INDEX_LOCAL(iTile,jTile,var)];
////                        //vars[var] = prim[j][i][var];
////                    }
////
////                    addSources(dU_dt,
////                               ucon, ucov, bcon, bcov,
////                               gcon, gcov, mhd, vars, 
////                               X1, X2);
////
//////                    // Compute fluxes along X1
//////                    X1 = i_TO_X1_FACE(i); X2 = j_TO_X2_CENTER(j);
//////                    ReconstructX1(primTile, iTile-1, jTile, primEdge, RIGHT);
//////                    ComputeFluxAndU(fluxR, uR, 
//////                                   ucon, ucov, bcon, bcov, 
//////                                   gcon, gcov, mhd, primEdge, 
//////                                   X1, X2, 1);
//////
//////                    ReconstructX1(primTile, iTile, jTile, primEdge, LEFT);
//////                    ComputeFluxAndU(fluxL, uL, 
//////                                    ucon, ucov, bcon, bcov, 
//////                                    gcon, gcov, mhd, primEdge, 
//////                                    X1, X2, 1);
//////
//////                    RiemannSolver(fluxR, fluxL, uR, uL, fluxX1L);
//////
//////
//////                    X1 = i_TO_X1_FACE(i+1); X2 = j_TO_X2_CENTER(j);
//////                    ReconstructX1(primTile, iTile, jTile, primEdge, RIGHT);
//////                    ComputeFluxAndU(fluxR, uR, 
//////                                    ucon, ucov, bcon, bcov, 
//////                                    gcon, gcov, mhd, primEdge, 
//////                                    X1, X2, 1);
//////
//////                    ReconstructX1(primTile, iTile+1, jTile, primEdge, LEFT);
//////                    ComputeFluxAndU(fluxL, uL, 
//////                                    ucon, ucov, bcon, bcov, 
//////                                    gcon, gcov, mhd, primEdge, 
//////                                    X1, X2, 1);
//////
//////                    RiemannSolver(fluxR, fluxL, uR, uL, fluxX1R);
//////
//////
//////                    // Compute fluxes along X2
//////                    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j);
//////                    ReconstructX2(primTile, iTile, jTile-1, primEdge, RIGHT);
//////                    ComputeFluxAndU(fluxR, uR, 
//////                                    ucon, ucov, bcon, bcov, 
//////                                    gcon, gcov, mhd, primEdge, 
//////                                    X1, X2, 2);
//////
//////                    ReconstructX2(primTile, iTile, jTile, primEdge, LEFT);
//////                    ComputeFluxAndU(fluxL, uL, 
//////                                    ucon, ucov, bcon, bcov, 
//////                                    gcon, gcov, mhd, primEdge, 
//////                                    X1, X2, 2);
//////
//////                    RiemannSolver(fluxR, fluxL, uR, uL, fluxX2L);
//////
//////
//////                    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j+1);
//////                    ReconstructX2(primTile, iTile, jTile, primEdge, RIGHT);
//////                    ComputeFluxAndU(fluxR, uR, 
//////                                    ucon, ucov, bcon, bcov, 
//////                                    gcon, gcov, mhd, primEdge, 
//////                                    X1, X2, 2);
//////
//////                    ReconstructX2(primTile, iTile, jTile+1, primEdge, LEFT);
//////                    ComputeFluxAndU(fluxL, uL, 
//////                                    ucon, ucov, bcon, bcov, 
//////                                    gcon, gcov, mhd, primEdge, 
//////                                    X1, X2, 2);
//////
//////                    RiemannSolver(fluxR, fluxL, uR, uL, fluxX2R);
////
////                    for (int var=0; var<DOF; var++) {
//////                        f[j][i][var] =  dU_dt[var] +
//////                                        (fluxX1R[var] - fluxX1L[var])/DX1 +
//////                                        (fluxX2R[var] - fluxX2L[var])/DX2;
////                        f[j][i][var] = dU_dt[var];
////                    if (isnan(f[j][i][var])) {
////                        printf("i = %d, j = %d, iTile = %d, jTile = %d, var=%d, f = %f\n",
////                                i, j, iTile, jTile, var, f[j][i][var]);
////                        printf("i = %d, j = %d, iTile = %d, jTile = %d, var=%d, prim= %f\n",
////                                i, j, iTile, jTile, var, prim[j][i][var]);
////                    }
////                    }
////
////                }
////            }
////        }
////    }
//
//    DMDAVecRestoreArrayDOF(dmda, PrimLocal, &prim);
//    DMDAVecRestoreArrayDOF(dmda, dPrim_dt, &dprim_dt);
//    DMDAVecRestoreArrayDOF(dmda, F, &f);
//                    
//    DMRestoreLocalVector(dmda, &PrimLocal);

    return(0.);
}


void Benchmark(TS ts, Vec Prim)
{
    PetscInt NIter = 1;
    std::clock_t start, end;
    PetscScalar time;

    PetscScalar t = 0.;
    Vec dPrim_dt, F;
    VecDuplicate(Prim, &dPrim_dt);
    VecDuplicate(Prim, &F);
    VecSet(dPrim_dt, 0.);
    VecSet(F, 0.);

    start = std::clock();
    for (int n=0; n < NIter; n++) {
        ComputeResidual(ts, t, Prim, dPrim_dt, F, NULL);
    }
    end = std::clock();

    time = (end - start)/(PetscScalar)CLOCKS_PER_SEC;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken for %d iterations = %f\n", NIter, time);

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD,"InitialResidual.h5",
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) F, "F");
    VecView(F, viewer);
    PetscViewerDestroy(&viewer);

    VecDestroy(&dPrim_dt);
    VecDestroy(&F);
}

void InitialConditionTest(TS ts, Vec X)
{
    DM dmda;
    TSGetDM(ts, &dmda);
    PetscScalar ***x;

    DMDAVecGetArrayDOF(dmda, X, &x);

    for (PetscInt j=0; j<N2; j++) {
        for (PetscInt i=0; i<N1; i++) {
            PetscScalar xcoord = DX1/2. + i*DX1;
            PetscScalar ycoord = DX2/2. + j*DX2;

            PetscScalar xcenter = 0.5;
            PetscScalar ycenter = 0.5;

            PetscScalar r = PetscSqrtScalar((xcoord-xcenter)*(xcoord-xcenter)+
                                            (ycoord-ycenter)*(ycoord-ycenter));
            
//            x[j][i][RHO] = 1. + exp(-r*r/0.01);
//
//            x[j][i][UU] = 1./(ADIABATIC_INDEX-1);
//            x[j][i][U1] = 4.95;
//            x[j][i][U2] = 4.95;
//            x[j][i][U3] = 0.;
//            x[j][i][B1] = 0.;
//            x[j][i][B2] = 0.;
//            x[j][i][B3] = 0.;

            x[j][i][RHO] = 25./(36.*M_PI);
            x[j][i][UU] = 5./(12.*M_PI*(ADIABATIC_INDEX - 1.));
            PetscScalar vx = -0.5*PetscSinScalar(2*M_PI*ycoord);
            PetscScalar vy = 0.5*PetscSinScalar(2*M_PI*xcoord);
            PetscScalar vz = 0.;
            PetscScalar gamma = 1./PetscSqrtScalar(1 - vx*vx - vy*vy - vz*vz);
            x[j][i][U1] = gamma*vx;
            x[j][i][U2] = gamma*vy;
            x[j][i][U3] = gamma*vz;
            x[j][i][B1] =
                -1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(2*M_PI*ycoord);
            x[j][i][B2] =
                1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(4*M_PI*xcoord);
            x[j][i][B3] = 0.0;

//                  Komissarov strong cylindrical explosion
//              PetscScalar R = 1.;
//              PetscScalar rho_outside = 1e-4, rho_inside = 1e-2;
//              PetscScalar u_outside = 3e-5/(5./3 - 1); 
//              PetscScalar u_inside = 1./(5./3 - 1);
//              PetscScalar alpha = 20.;
//              PetscScalar norm_factor = (1. - tanh(-R))/2.;
//
//              x[j][i][RHO] = ((rho_inside - rho_outside)*
//                              (1. - tanh(pow(r, alpha)-R))
//                              /2./norm_factor + rho_outside);
//
//              x[j][i][UU] = ((u_inside - u_outside)*
//                            (1. - tanh(pow(r, alpha)-R))
//                            /2./norm_factor + u_outside);
//
//              PetscScalar vx = 0.;
//              PetscScalar vy = 0.;
//              PetscScalar vz = 0.;
//              PetscScalar gamma = 1.;
//              x[j][i][U1] = gamma*vx;
//              x[j][i][U2] = gamma*vy;
//              x[j][i][U3] = gamma*vz;
//              x[j][i][B1] = .1;
//              x[j][i][B2] = 0.;
//              x[j][i][B3] = 0.0;

        }
    }

    DMDAVecRestoreArrayDOF(dmda, X, &x);
}

void InitialCondition(TS ts, Vec Prim)
{
    DM dmda;
    PetscScalar ***prim, r, theta;
    PetscScalar X1, X2;

    PetscScalar l, delta, sigma, A, lnOfh;
    PetscScalar thetaIn, deltaIn, sigmaIn, AIn;
    PetscScalar uPhiTmp, uconBL[NDIM], uconKS[NDIM], uconKSPrime[NDIM]; 
    PetscScalar gcovBL[NDIM][NDIM], gconBL[NDIM][NDIM], transform[NDIM][NDIM];
    PetscScalar rhoMax = 0., uMax = 0., bsqrMax = 0., rhoAv=0.;
    PetscScalar q, norm, betaActual;
    PetscScalar AA, BB, CC, DD, discriminent, mu;

    int X1Start, X2Start;
    int X1Size, X2Size;

    PetscScalar randNum;
    PetscRandom randNumGen;
    PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
    PetscRandomSetType(randNumGen, PETSCRAND48);

    TSGetDM(ts, &dmda);

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    PetscScalar AVector[X2Size+2*NG][X1Size+2*NG];

    Vec localPrim;
    DMGetLocalVector(dmda, &localPrim);
    DMDAVecGetArrayDOF(dmda, localPrim, &prim);

	l =  ( ( (pow(A_SPIN, 2) - 2.*A_SPIN*sqrt(R_MAX) + pow(R_MAX, 2.)) *
		     ( (-2.*A_SPIN*R_MAX*(pow(A_SPIN, 2.) - 2.*A_SPIN*sqrt(R_MAX) +
		        pow(R_MAX, 2.))) / sqrt(2.*A_SPIN*sqrt(R_MAX) + (-3. +
                R_MAX)*R_MAX) +
		       ((A_SPIN + (-2. + R_MAX)*sqrt(R_MAX))*(pow(R_MAX, 3) +
				pow(A_SPIN, 2)*(2. + R_MAX))) / sqrt(1 + (2.*A_SPIN) /
			   pow(R_MAX, 1.5) - 3./R_MAX) ) ) / \
           (pow(R_MAX, 3)*sqrt(2.*A_SPIN*sqrt(R_MAX) +
           (-3.+R_MAX)*R_MAX)*(pow(A_SPIN,2) + (-2.+R_MAX)*R_MAX)) );

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++)
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {
            
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            BLCoords(&r, &theta, X1, X2);

            delta = r*r - 2*M*r + A_SPIN*A_SPIN;
            sigma = r*r + A_SPIN*A_SPIN*cos(theta)*cos(theta);
            A = (r*r + A_SPIN*A_SPIN)*(r*r + A_SPIN*A_SPIN) - 
                delta*A_SPIN*A_SPIN*sin(theta)*sin(theta);

            thetaIn = M_PI/2.;

            deltaIn = R_MIN*R_MIN - 2*M*R_MIN + A_SPIN*A_SPIN;
            sigmaIn = R_MIN*R_MIN + A_SPIN*A_SPIN*cos(thetaIn)*cos(thetaIn);
            AIn = (R_MIN*R_MIN + A_SPIN*A_SPIN)*(R_MIN*R_MIN + A_SPIN*A_SPIN) - 
                  deltaIn*A_SPIN*A_SPIN*sin(thetaIn)*sin(thetaIn);

            if (r >=R_MIN) {

			    lnOfh = 0.5*log((1. + sqrt(1. + 4.*(l*l*sigma*sigma)*delta/\
                                (A*sin(theta)*A*sin(theta))))/(sigma*delta/A))-
			            0.5*sqrt(1. + 4.*(l*l*sigma*sigma)*delta /
					             (A*A*sin(theta)*sin(theta))) -2.*A_SPIN*r*l/A-
			           (0.5*log((1. + sqrt(1. +
				                           4.*(l*l*sigmaIn*sigmaIn)*deltaIn/ \
				        (AIn*AIn*sin(thetaIn)*sin(thetaIn)))) /
				        (sigmaIn * deltaIn / AIn)) - 0.5 * sqrt(1. +
					    4.*(l*l*sigmaIn*sigmaIn)*deltaIn/ \
                        (AIn*AIn*sin(thetaIn)*sin(thetaIn))) - 
                        2.*A_SPIN*R_MIN*l/AIn);
		    } else {
			    lnOfh = 1.;
            }

            if (lnOfh <0. || r < R_MIN) {

                prim[j][i][RHO] = 1e-7*RHO_MIN;
                prim[j][i][UU] = 1e-7*U_MIN;
                prim[j][i][U1] = 0.;
                prim[j][i][U2] = 0.;
                prim[j][i][U3] = 0.;

            } else {

                prim[j][i][RHO] = pow((exp(lnOfh) - 1.)*
                    (ADIABATIC_INDEX-1.)/(KAPPA*ADIABATIC_INDEX),
                    1./(ADIABATIC_INDEX-1.));

                if (prim[j][i][RHO] > rhoMax)
                    rhoMax = prim[j][i][RHO];

                PetscRandomGetValue(randNumGen, &randNum);
                prim[j][i][UU] = KAPPA*pow(prim[j][i][RHO], ADIABATIC_INDEX)/\
                              (ADIABATIC_INDEX-1.)*(1. + 4e-2*(randNum-0.5));

                if (prim[j][i][UU] > uMax && r > R_MIN)
                    uMax = prim[j][i][UU];

                uPhiTmp = sqrt((-1. + sqrt(1. + 4.*l*l*(sigma*sigma*delta/\
                               (A*A*sin(theta)*sin(theta))))) / 2.);
			
                uconBL[1] = 0.;
                uconBL[2] = 0.;
                uconBL[3] = 2.*A_SPIN*r*sqrt(1. + uPhiTmp*uPhiTmp)/ \
                            sqrt(A*sigma*delta) + 
                            sqrt(sigma/A)*uPhiTmp/sin(theta);
                
                // Transform uconBoyerLinquist to uconKerrSchild and set to prim

                for (int ii=0; ii<NDIM; ii++)
                    for (int jj=0; jj<NDIM; jj++) {
                        gcovBL[ii][jj] = 0.;
                        gconBL[ii][jj] = 0.;
                        transform[ii][jj] = 0.;
                    }

                DD = 1. - 2./r + A_SPIN*A_SPIN/(r*r);
                mu = 1 + A_SPIN*A_SPIN*cos(theta)*cos(theta)/(r*r);

                gcovBL[0][0] = -(1. - 2./(r*mu));
                gcovBL[0][3] = -2.*A_SPIN*sin(theta)*sin(theta)/(r*mu);
                gcovBL[3][0] = gcovBL[0][3];
                gcovBL[1][1] = mu/DD;
                gcovBL[2][2] = r*r*mu;
                gcovBL[3][3] = r*r*sin(theta)*sin(theta)*\
                               (1. + A_SPIN*A_SPIN/(r*r) + 
                                2.*A_SPIN*A_SPIN*sin(theta)*sin(theta)/\
                                (r*r*r*mu));

                gconBL[0][0] = -1. -2.*(1 + A_SPIN*A_SPIN/(r*r))/(r*DD*mu);
                gconBL[0][3] = -2.*A_SPIN/(r*r*r*DD*mu);
                gconBL[3][0] = gconBL[0][3];
                gconBL[1][1] = DD/mu;
                gconBL[2][2] = 1./(r*r*mu);
                gconBL[3][3] = (1. - 2./(r*mu))/(r*r*sin(theta)*sin(theta)*DD);

                transform[0][0] = 1.;
                transform[1][1] = 1.;
                transform[2][2] = 1.;
                transform[3][3] = 1.;
                transform[0][1] = 2.*r/(r*r* - 2.*r + A_SPIN*A_SPIN);
                transform[3][1] = A_SPIN/(r*r - 2.*r + A_SPIN*A_SPIN);

                AA = gcovBL[0][0];
                BB = 2.*(gcovBL[0][1]*uconBL[1] + 
                         gcovBL[0][2]*uconBL[2] +
                         gcovBL[0][3]*uconBL[3]);
                CC = 1. + gcovBL[1][1]*uconBL[1]*uconBL[1] +
                          gcovBL[2][2]*uconBL[2]*uconBL[2] +
                          gcovBL[3][3]*uconBL[3]*uconBL[3] +
                      2.*(gcovBL[1][2]*uconBL[1]*uconBL[2] +
                          gcovBL[1][3]*uconBL[1]*uconBL[3] +
                          gcovBL[2][3]*uconBL[2]*uconBL[3]);

                discriminent = BB*BB - 4.*AA*CC;
                uconBL[0] = -(BB + sqrt(discriminent))/(2.*AA);

                uconKS[0] = transform[0][0]*uconBL[0] + 
                            transform[0][1]*uconBL[1] +
                            transform[0][2]*uconBL[2] +
                            transform[0][3]*uconBL[3];

                uconKS[1] = transform[1][0]*uconBL[0] + 
                            transform[1][1]*uconBL[1] +
                            transform[1][2]*uconBL[2] +
                            transform[1][3]*uconBL[3];

                uconKS[2] = transform[2][0]*uconBL[0] + 
                            transform[2][1]*uconBL[1] +
                            transform[2][2]*uconBL[2] +
                            transform[2][3]*uconBL[3];

                uconKS[3] = transform[3][0]*uconBL[0] + 
                            transform[3][1]*uconBL[1] +
                            transform[3][2]*uconBL[2] +
                            transform[3][3]*uconBL[3];

                PetscScalar rFactor, hFactor;
                rFactor = r - R0;
                hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*X2);
                uconKSPrime[0] = uconKS[0];
                uconKSPrime[1] = uconKS[1]*(1./rFactor);
                uconKSPrime[2] = uconKS[2]*(1./hFactor);
                uconKSPrime[3] = uconKS[3];

                PetscScalar gcov[NDIM][NDIM], gcon[NDIM][NDIM];
                PetscScalar gdet, alpha;
                gCovCalc(gcov, X1, X2);
                gDetCalc(&gdet, gcov);
                gConCalc(gcon, gcov, gdet);
                alphaCalc(&alpha, gcon);

                prim[j][i][U1] = uconKSPrime[1] + 
                                 alpha*alpha*uconKSPrime[0]*gcon[0][1];
                prim[j][i][U2] = uconKSPrime[2] + 
                                 alpha*alpha*uconKSPrime[0]*gcon[0][2];
                prim[j][i][U3] = uconKSPrime[3] +
                                 alpha*alpha*uconKSPrime[0]*gcon[0][3];

            }

            prim[j][i][B1] = 0.;
            prim[j][i][B2] = 0.;
            prim[j][i][B3] = 0.;

        }


    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        
            prim[j][i][RHO] = prim[j][i][RHO]/rhoMax;
            prim[j][i][UU] = prim[j][i][UU]/rhoMax;

            AVector[j][i] = 0.;
        }
    }

    uMax = uMax/rhoMax;

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            REAL gcov[NDIM][NDIM];
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloor(vars, gamma, X1, X2);

            for (int var=0; var<DOF; var++)
                prim[j][i][var] = vars[var];
        }
    }

    for (int j=X2Start-1; j<X2Start+X2Size+1; j++) {
        for (int i=X1Start-1; i<X1Start+X1Size+1; i++) {        
            rhoAv = 0.25*(prim[j][i][RHO] + prim[j][i-1][RHO] +
                          prim[j-1][i][RHO] + prim[j-1][i-1][RHO]);
            
            q = rhoAv - 0.2;

            if (q > 0.)
                AVector[j][i] = q;

        }
    }

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
            
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            PetscScalar gcov[NDIM][NDIM], gcon[NDIM][NDIM];
            PetscScalar g, gdet, alpha;
            gCovCalc(gcov, X1, X2);
            gDetCalc(&gdet, gcov);
            gConCalc(gcon, gcov, gdet);
            alphaCalc(&alpha, gcon);
            g = sqrt(-gdet);

            prim[j][i][B1] = -(AVector[j][i] - AVector[j+1][i] +
                               AVector[j][i+1] - AVector[j+1][i+1])/\
                              (2.*DX2*g);

            prim[j][i][B2] = (AVector[j][i] + AVector[j+1][i] -
                            AVector[j][i+1] - AVector[j+1][i+1])/\
                            (2.*DX1*g);

            prim[j][i][B3] = 0.;

            PetscScalar gamma, ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
            PetscScalar bsqr, var[DOF];

            for (int n=0; n<DOF; n++)
                var[n] = prim[j][i][n];

            gammaCalc(&gamma, var, gcov);
            uconCalc(ucon, gamma, alpha, var, gcon);
            covFromCon(ucov, ucon, gcov);
            bconCalc(bcon, var, ucon, ucov);
            covFromCon(bcov, bcon, gcov);
            bSqrCalc(&bsqr, bcon, bcov);

            if (bsqr > bsqrMax)
                bsqrMax = bsqr;

        }
    }
    betaActual = (ADIABATIC_INDEX-1.)*uMax/(0.5*bsqrMax);
    norm = sqrt(betaActual/BETA);

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
            prim[j][i][B1] = prim[j][i][B1]*norm;
            prim[j][i][B2] = prim[j][i][B2]*norm;
        }
    }

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            REAL gcov[NDIM][NDIM];
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloor(vars, gamma, X1, X2);

            prim[j][i][RHO] = (vars[RHO]);
            prim[j][i][UU] = (vars[UU]);

            for (int var=2; var<DOF; var++)
                prim[j][i][var] = vars[var];

        }
    }

    DMLocalToGlobalBegin(dmda, localPrim, INSERT_VALUES, Prim);
    DMLocalToGlobalEnd(dmda, localPrim, INSERT_VALUES, Prim);

    DMDAVecRestoreArrayDOF(dmda, localPrim, &prim);
    DMRestoreLocalVector(dmda, &localPrim);

    PetscRandomDestroy(&randNumGen);
}


PetscErrorCode Monitor(TS ts, 
                       PetscInt step,
                       PetscReal t,
                       Vec Prim,
                       void *ptr)
{
    DM dmda;
    TSGetDM(ts, &dmda);
    PetscScalar ***prim;
    PetscScalar gcov[NDIM][NDIM];

    DMDAVecGetArrayDOF(dmda, Prim, &prim);

    for (int j=0; j<N2; j++) {
        for (int i=0; i<N1; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloor(vars, gamma, X1, X2);

            for (int var=0; var<DOF; var++)
                prim[j][i][var] = vars[var];
        }
    }

    DMDAVecRestoreArrayDOF(dmda, Prim, &prim);

    REAL dt, dtDump;
    TSGetTimeStep(ts, &dt);
    
    static PetscInt counter = 1;
    static PetscScalar tDump = 0.;
    dtDump = 1.;

    if (t > tDump) {
        printf("Dumping data..\n");
        std::string filename = "plot";
        filename += std::to_string(counter);
        filename += ".h5";

        PetscViewer viewer;
        PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(),
                            FILE_MODE_WRITE, &viewer);
        PetscObjectSetName((PetscObject) Prim, "soln");
        VecView(Prim, viewer);
        PetscViewerDestroy(&viewer);

        tDump = tDump + dtDump;
        counter++;
    }


//    if (t>250) {
//        dt = 0.01;
//        TSSetTimeStep(ts, dt);
//    }
//    if (dt<0.025) {
//        dt = 0.025;
//        TSSetTimeStep(ts, dt);
//    }
//    TSSetTimeStep(ts, 0.01);


    return(0.);
}
