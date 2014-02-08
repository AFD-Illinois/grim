#include "grim.h"

void BLCoord(REAL X1, REAL X2, REAL r, REAL theta) 
{
    r = exp(X1) + R0;
    theta = M_PI*X2 + 
            ((1. - H_SLOPE)/2.)*sin(2.*M_PI*X2);

}

PetscErrorCode SetCoordinates(DM dmda,
                              PetscInt X1Start, PetscInt X1Size,
                              PetscInt X2Start, PetscInt X2Size)
{
    DM coordDM;
    Vec coordVec;

    DMDASetUniformCoordinates(dmda,0.0,1.0,0.0,1.0,0.0,1.0);
    DMGetCoordinateDM(dmda, &coordDM);
    DMGetCoordinatesLocal(dmda, &coordVec);

    PetscScalar r, theta;
    DMDACoor2d **coord;
    
    DMDAVecGetArray(coordDM, coordVec, &coord);
//
//    for (int j=X2Start; j<X2Start+X2Size; j++) {
//        for (int i=X1Start; i<X1Start+X1Size; i++) {
//
//            PetscScalar X1 = X1_MIN + DX1*(i+X1Start);
//            PetscScalar X2 = X2_MIN + DX2*(j+X2Start);
//
//            BLCoord(X1, X2, r, theta);
//
//            coord[j][i].x = r*PetscSinScalar(theta);
//            coord[j][i].y = r*PetscCosScalar(theta);
//        }
//    }
    DMDAVecRestoreArray(coordDM, coordVec, &coord);

    DMSetCoordinatesLocal(dmda, coordVec);
}

int main(int argc, char **argv)
{
    TS ts;
    Vec soln, coordVec;
    DM da, coordDM;
    DM DAs[DOF];

    int X1Start, X2Start;
    int X1Size, X2Size;

    PetscInitialize(&argc, &argv, PETSC_NULL, help);

//    viennacl::ocl::set_context_platform_index(0, 1);

    DMCompositeCreate(PETSC_COMM_WORLD, &da);

    for (int var=0; var < DOF; var++) {
        DMDACreate2d(PETSC_COMM_WORLD, 
                     DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                     DMDA_STENCIL_STAR,
                     N1, N2,
                     PETSC_DECIDE, PETSC_DECIDE,
                     1, NG, PETSC_NULL, PETSC_NULL, &DAs[var]);

        DMDAGetCorners(DAs[var], 
                       &X1Start, &X2Start, NULL,
                       &X1Size, &X2Size, NULL);

        SetCoordinates(DAs[var],
                       X1Start, X1Size, X2Start, X2Size);

        DMCompositeAddDM(da, DAs[var]);
       
        DMDestroy(&DAs[var]);
    }

    DMCreateGlobalVector(da, &soln);

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetDM(ts, da);
    TSSetIFunction(ts, PETSC_NULL, ComputeResidual, NULL);

    std::ifstream sourceFile("computeresidual.cl");
    std::string sourceCode((std::istreambuf_iterator<char>(sourceFile)),
                            std::istreambuf_iterator<char>());

    std::string BuildOptions("\
                              -D X1_SIZE=" +
                             std::to_string(X1Size) +
                             " -D X2_SIZE=" + 
                             std::to_string(X2Size) +
                             " -D TOTAL_X1_SIZE=" + 
                             std::to_string(X1Size+2*NG) + 
                             " -D TOTAL_X2_SIZE=" +
                             std::to_string(X2Size+2*NG));

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

//    kernel.global_work_size(N2, N1);
//    kernel.local_work_size(TILE_SIZE_X2, TILE_SIZE_X1);

//    clErr = cl::Platform::get(&platforms);
//    CheckCLErrors(clErr, "cl::Platform::get");
//
//    clErr = platforms.at(1).getDevices(CL_DEVICE_TYPE_CPU, &devices);
//    CheckCLErrors(clErr, "cl::Platform::getDevices");
//
//    context = cl::Context(devices, NULL, NULL, NULL, &clErr);
//    CheckCLErrors(clErr, "cl::Context::Context");
//
//    queue = cl::CommandQueue(context, devices.at(0), 0, &clErr);
//    CheckCLErrors(clErr, "cl::CommandQueue::CommandQueue");
//
//    std::ifstream sourceFile("computeresidual.cl");
//    std::string sourceCode((std::istreambuf_iterator<char>(sourceFile)),
//                            std::istreambuf_iterator<char>());
//    cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(),
//                                sourceCode.length()+1));
//    
//    program = cl::Program(context, source, &clErr);
//    CheckCLErrors(clErr, "cl::Program::Program");
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
//    PetscScalar start = std::clock();
//    clErr = program.build(devices, BuildOptions.c_str(), NULL, NULL);
//    const char *buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(
//                                                devices.at(0),
//                                                &clErr).c_str();
//    PetscPrintf(PETSC_COMM_WORLD, "%s\n", buildlog);
//    CheckCLErrors(clErr, "cl::Program::build");
//    PetscScalar end = std::clock();
//
//    PetscScalar time = (end - start)/(PetscScalar)CLOCKS_PER_SEC;
//    PetscPrintf(PETSC_COMM_WORLD, 
//                "Time taken for kernel compilation = %f\n", time);
//
//
//    kernel = cl::Kernel(program, "ComputeResidual", &clErr);
//    CheckCLErrors(clErr, "cl::Kernel::Kernel");

    InitialCondition(ts, soln);

//    Benchmark(ts, soln);

    TSSetSolution(ts, soln);
    TSSetType(ts, TSTHETA);
    TSSetFromOptions(ts);

    TSSolve(ts, soln);

//    VecView(soln, PETSC_VIEWER_STDOUT_WORLD);
//    PetscViewer viewer;
//    PetscViewerHDF5Open(PETSC_COMM_WORLD,"soln.h5", FILE_MODE_WRITE, &viewer);
//    PetscObjectSetName((PetscObject) soln, "var");
//    VecView(soln, viewer);
//    PetscViewerDestroy(&viewer);    


    DMDestroy(&da);
    VecDestroy(&soln);
    TSDestroy(&ts);

    PetscFinalize();
    return(0);
}

int counter = 0;
PetscErrorCode ComputeResidual(TS ts,
                               PetscScalar t,
                               Vec Prim, Vec dPrim_dt,
                               Vec F, void *ptr)
{
//    counter++;
////    printf("In ComputeResidual %d\n", counter);
//    const viennacl::vector<PetscScalar> *vars, *dvars_dt;
//    viennacl::vector<PetscScalar> *Fvars;
//
//    VecViennaCLGetArrayRead(Prim, &vars);
//    VecViennaCLGetArrayRead(dPrim_dt, &dvars_dt);
//    VecViennaCLGetArrayWrite(F, &Fvars);
//
//    viennacl::ocl::kernel kernel = program.get_kernel("ComputeResidual");
//
//    kernel.local_work_size(0, TILE_SIZE_X1);
//    kernel.local_work_size(1, TILE_SIZE_X2);
//    kernel.global_work_size(0, N1);
//    kernel.global_work_size(1, N2);
//
////    viennacl::ocl::enqueue(kernel(*vars, *dvars_dt, *Fvars));
////    viennacl::backend::finish();
//
//    VecViennaCLRestoreArrayRead(Prim, &vars);
//    VecViennaCLRestoreArrayRead(dPrim_dt, &dvars_dt);
//    VecViennaCLRestoreArrayWrite(F, &Fvars);

//    PetscScalar *prim, *dprim_dt, *f;
//    VecGetArray(Prim, &prim);
//    VecGetArray(dPrim_dt, &dprim_dt);
//    VecGetArray(F, &f);
//    
//    for (int var=0; var<DOF; var++)
//        for (int j=0; j<N2-1; j++)
//            for (int i=0; i<N1-1; i++) {
//                f[INDEX_GLOBAL(i,j,var)] = dprim_dt[INDEX_GLOBAL(i,j,var)] -
//                                         (prim[INDEX_GLOBAL(i+1,j,var)] - 
//                                          prim[INDEX_GLOBAL(i,j,var)])/DX1 -
//                                         (prim[INDEX_GLOBAL(i,j+1,var)] - 
//                                          prim[INDEX_GLOBAL(i,j,var)])/DX2;
//            }
//    VecRestoreArray(Prim, &prim);
//    VecRestoreArray(dPrim_dt, &dprim_dt);
//    VecRestoreArray(F, &f);

//    DM da;
//    DM DAs[DOF];
//    Vec Vecs[DOF], dVecs_dt[DOF]; //Local vectors with ghost zones.
//    Vec Fvecs[DOF]; //Local vectors without ghost zones.
//    PetscScalar **vars[DOF], **dvars_dt[DOF], **Fvars[DOF];
//    int X1Start, X2Start;
//    int X1Size, X2Size;
//
//    TSGetDM(ts, &da);
//
//    DMCompositeGetEntriesArray(da, DAs);
//    DMDAGetCorners(DAs[0], 
//                   &X1Start, &X2Start, NULL,
//                   &X1Size, &X2Size, NULL);
//
//    for (int var=0; var<DOF; var++) {
//        DMGetLocalVector(DAs[var], &Vecs[var]);
//        DMGetLocalVector(DAs[var], &dVecs_dt[var]);
//    }
//    
//    // DMCompositeScatterArray takes A LOT OF TIME! Find out why? Possible
//    // performance bug.
//    DMCompositeScatterArray(da, Prim, Vecs);
//    DMCompositeScatterArray(da, dPrim_dt, dVecs_dt);
//    DMCompositeGetAccessArray(da, F, DOF, NULL, Fvecs);
//
//    for (int var=0; var<DOF; var++) {
//        DMDAVecGetArray(DAs[var], Vecs[var], &vars[var]);
//        DMDAVecGetArray(DAs[var], dVecs_dt[var], &dvars_dt[var]);
//        DMDAVecGetArray(DAs[var], Fvecs[var], &Fvars[var]);
//    }
//
//    for (int var=0; var<DOF; var++)
//    for (int j=X2Start; j<X2Start+X2Size; j++)
//        for (int i=X1Start; i<X1Start+X1Size; i++) {
//            Fvars[var][j][i] = dvars_dt[var][j][i] +
//                               (vars[var][j][i+1] - vars[var][j][i])/DX1 +
//                               (vars[var][j+1][i] - vars[var][j][i])/DX2;
//        }
//
//    for (int var=0; var<DOF; var++) {
//        DMDAVecRestoreArray(DAs[var], Vecs[var], &vars[var]);
//        DMDAVecRestoreArray(DAs[var], dVecs_dt[var], &dvars_dt[var]);
//        DMDAVecRestoreArray(DAs[var], Fvecs[var], &Fvars[var]);
//
//        DMRestoreLocalVector(DAs[var], &Vecs[var]);
//        DMRestoreLocalVector(DAs[var], &dVecs_dt[var]);
//    }
//
//    DMCompositeRestoreAccessArray(da, F, DOF, NULL, Fvecs);



    return(0.);
}

//PetscErrorCode ComputeResidual(TS ts,
//                               PetscScalar t,
//                               Vec Prim, Vec dPrim_dt,
//                               Vec F, void *ptr)
//{
//    DM da;
//    DM DAs[DOF];
//    Vec Vecs[DOF], dVecs_dt[DOF]; //Local vectors with ghost zones.
//    Vec Fvecs[DOF]; //Local vectors without ghost zones.
//    PetscScalar **vars[DOF], **dvars_dt[DOF], **Fvars[DOF];
//
//    int X1Start, X2Start;
//    int X1Size, X2Size;
//    int totalSize, sizeWithoutNG, totalTileSize;
//
//    TSGetDM(ts, &da);
//
//    DMCompositeGetEntriesArray(da, DAs);
//    DMDAGetCorners(DAs[0], 
//                   &X1Start, &X2Start, NULL,
//                   &X1Size, &X2Size, NULL);
//
//    totalSize = (X1Size + 2*NG)*(X2Size + 2*NG)*sizeof(PetscScalar);
//    sizeWithoutNG = X1Size*X2Size*sizeof(PetscScalar);
//    totalTileSize = sizeof(PetscScalar)*(TILE_SIZE+2*NG)*DOF*(2*NG+1);
//    
//    for (int var=0; var<DOF; var++) {
//        DMGetLocalVector(DAs[var], &Vecs[var]);
//        DMGetLocalVector(DAs[var], &dVecs_dt[var]);
//    }
//    
//    // DMCompositeScatterArray takes A LOT OF TIME! Find out why? Possible
//    // performance bug.
//    DMCompositeScatterArray(da, Prim, Vecs);
//    DMCompositeScatterArray(da, dPrim_dt, dVecs_dt);
//    DMCompositeGetAccessArray(da, F, DOF, NULL, Fvecs);
//
//    for (int var=0; var<DOF; var++) {
//        DMDAVecGetArray(DAs[var], Vecs[var], &vars[var]);
//        DMDAVecGetArray(DAs[var], dVecs_dt[var], &dvars_dt[var]);
//        DMDAVecGetArray(DAs[var], Fvecs[var], &Fvars[var]);
//    }
//
//    cl::Buffer PrimBuffers[DOF], dPrimBuffers_dt[DOF], Fbuffers[DOF];
//
//    for (int var=0; var<DOF; var++) {
//        PrimBuffers[var] = cl::Buffer(context,
//                                      CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
//                                      totalSize, 
//                                      &(vars[var][X2Start-NG][X1Start-NG]),
//                                      &clErr);
//        dPrimBuffers_dt[var] = cl::Buffer(context,
//                                          CL_MEM_USE_HOST_PTR | 
//                                          CL_MEM_READ_ONLY,
//                                          totalSize, 
//                                          &(dvars_dt[var][X2Start-NG][X1Start-NG]),
//                                          &clErr);
//        Fbuffers[var] = cl::Buffer(context,
//                                   CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
//                                   sizeWithoutNG, 
//                                   &(Fvars[var][X2Start][X1Start]),
//                                   &clErr);
//
//        clErr = kernel.setArg(var, PrimBuffers[var]);
//        clErr = kernel.setArg(var+DOF, dPrimBuffers_dt[var]);
//        clErr = kernel.setArg(var+2*DOF, Fbuffers[var]);
//    }
//
//    // WARNING! Local memory MUST NOT exceed CL_DEVICE_LOCAL_MEM_SIZE. Else code
//    // will fail silently. Put a check for this in main().
//    cl::LocalSpaceArg primTile = cl::__local(totalTileSize);
//
//    clErr = kernel.setArg(3*DOF, primTile);
//
//    cl::NDRange global(X1Size);
//    cl::NDRange local(TILE_SIZE);
//    clErr = queue.enqueueNDRangeKernel(kernel,
//                                       cl::NullRange,
//                                       global, local,
//                                       NULL, NULL);
//
//    for (int var=0; var<DOF; var++) {
//        Fvars[var][X2Start] = 
//                    (PetscScalar*)queue.enqueueMapBuffer(Fbuffers[var],
//                                                         CL_FALSE,
//                                                         CL_MAP_READ,
//                                                         0, sizeWithoutNG,
//                                                         NULL, NULL, &clErr);
//    }
//
//    clErr = queue.finish();
//
//    for (int var=0; var<DOF; var++) {
//        DMDAVecRestoreArray(DAs[var], Vecs[var], &vars[var]);
//        DMDAVecRestoreArray(DAs[var], dVecs_dt[var], &dvars_dt[var]);
//        DMDAVecRestoreArray(DAs[var], Fvecs[var], &Fvars[var]);
//
//        DMRestoreLocalVector(DAs[var], &Vecs[var]);
//        DMRestoreLocalVector(DAs[var], &dVecs_dt[var]);
//    }
//
//    DMCompositeRestoreAccessArray(da, F, DOF, NULL, Fvecs);
//
//    return(0.);
//}


void Benchmark(TS ts, Vec Prim)
{
    PetscInt NIter = 1;
    std::clock_t start, end;
    PetscScalar time;

    PetscScalar t = 0.;
    Vec dPrim_dt, F;
    VecDuplicate(Prim, &dPrim_dt);
    VecDuplicate(Prim, &F);
    VecSet(dPrim_dt, 1.);
    VecSet(F, 2.);

    start = std::clock();
    for (int n=0; n < NIter; n++) {
        printf("n = %d\n", n);
        ComputeResidual(ts, t, Prim, dPrim_dt, F, NULL);
    }
    end = std::clock();

    time = (end - start)/(PetscScalar)CLOCKS_PER_SEC;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken for %d iterations = %f\n", NIter, time);

    VecDestroy(&dPrim_dt);
    VecDestroy(&F);
}

void InitialCondition(TS ts, Vec Prim)
{
    DM da;
    DM DAs[DOF];
    Vec Vecs[DOF];
    PetscScalar **vars[DOF];

    int X1Start, X2Start;
    int X1Size, X2Size;
    int X1End, X2End;

    TSGetDM(ts, &da);

    DMCompositeGetEntriesArray(da, DAs);
    for (int var=0; var<DOF; var++) {
        DMGetLocalVector(DAs[var], &Vecs[var]);
        DMDAVecGetArray(DAs[var], Vecs[var], &vars[var]);
    }

    DMDAGetCorners(DAs[0], 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    X1End = X1Start + X1Size;
    X2End = X2Start + X2Size;

    for (PetscInt j=X2Start; j<X2End; j++) {
        for (PetscInt i=X1Start; i<X1End; i++) {
            
            PetscScalar X1Coord = X1_MIN + DX1/2. + i*DX1;
            PetscScalar X2Coord = X2_MIN + DX2/2. + j*DX2;

            PetscScalar X1Center = (X1_MIN + X1_MAX)/2.;
            PetscScalar X2Center = (X2_MIN + X2_MAX)/2.;

            PetscScalar r = PetscSqrtScalar(
                                PetscPowScalar(X1Coord-X1Center, 2.0) + 
                                PetscPowScalar(X2Coord-X2Center, 2.0));

// Magnetized field loop advection
            PetscScalar A0 = 1e-3, R = 0.2;
            PetscScalar vx = 0.2/PetscSqrtScalar(6);
            PetscScalar vy = 0.1/PetscSqrtScalar(6);
            PetscScalar vz = 0.1/PetscSqrtScalar(6);
            PetscScalar gamma = 1./PetscSqrtScalar(1 - vx*vx - vy*vy - vz*vz);

//            vars[RHO][j][i] = 1.;
//            vars[UU][j][i] = 3./(5./3 - 1.);
//            vars[U1][j][i] = gamma*vx;
//            vars[U2][j][i] = gamma*vy;
//            vars[U3][j][i] = gamma*vz;
//            vars[B1][j][i] = -2*A0*X2Coord/(R*R)*PetscExpScalar(-r*r/(R*R));
//            vars[B2][j][i] = 2*A0*X1Coord/(R*R)*PetscExpScalar(-r*r/(R*R));
//            vars[B3][j][i] = 0.;

            for (int var=0; var<DOF; var++)
                vars[var][j][i] = exp(-r*r/.01);

//            vars[RHO][j][i] = 2.2;
//            vars[UU][j][i] = 2.1;
//            vars[U1][j][i] = 2.5;
//            vars[U2][j][i] = 2.3;
//            vars[U3][j][i] = 2.4;
//            vars[B1][j][i] = 2.9;
//            vars[B2][j][i] = 2.8;
//            vars[B3][j][i] = 2.7;
        }
    }

    DMCompositeGatherArray(da, Prim, INSERT_VALUES, Vecs);

    for (int var=0; var<DOF; var++) {
        DMDAVecRestoreArray(DAs[var], Vecs[var], &vars[var]);
        DMRestoreLocalVector(DAs[var], &Vecs[var]);
    }
}


//PetscErrorCode ComputeFlux2D(PetscScalar ***flux, PetscScalar ***x,
//                             PetscInt dir, AppCtx *ctx)
//{
//    PetscInt i, j;
//    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;
//
//    for (j=ctx->ystart-2; j<ctx->ystart + ctx->ysize+2; j++) 
//        for (i=ctx->xstart-2; i<ctx->xstart + ctx->xsize+2; i++) {
//            P = (5./3 - 1)*x[j][i][u];
//
//            gamma = PetscSqrtScalar(1. + x[j][i][u1]*x[j][i][u1] + x[j][i][u2]*x[j][i][u2] +
//                                         x[j][i][u3]*x[j][i][u3]);
//            u_con[0] = gamma;
//            u_con[1] = x[j][i][u1];
//            u_con[2] = x[j][i][u2];
//            u_con[3] = x[j][i][u3];
//
//            u_cov[0] = -gamma;
//            u_cov[1] = x[j][i][u1];
//            u_cov[2] = x[j][i][u2];
//            u_cov[3] = x[j][i][u3];
//
//            b_con[0] = x[j][i][B1]*u_con[1] + x[j][i][B2]*u_con[2] + x[j][i][B3]*u_con[3];
//            b_con[1] = (x[j][i][B1] + b_con[0]*u_con[1])/u_con[0];
//            b_con[2] = (x[j][i][B2] + b_con[0]*u_con[2])/u_con[0];
//            b_con[3] = (x[j][i][B3] + b_con[0]*u_con[3])/u_con[0];
//
//            b_cov[0] = -b_con[0];
//            b_cov[1] = b_con[1];
//            b_cov[2] = b_con[2];
//            b_cov[3] = b_con[3];
//
//            b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
//                     b_con[3]*b_con[3];
//
//            flux[j][i][rho] = x[j][i][rho]*u_con[dir];
//        
//            flux[j][i][B1] = b_con[1]*u_con[dir] - b_con[dir]*u_con[1];
//            flux[j][i][B2] = b_con[2]*u_con[dir] - b_con[dir]*u_con[2];
//            flux[j][i][B3] = b_con[3]*u_con[dir] - b_con[dir]*u_con[3];
//
//            flux[j][i][u] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[0] + 
//                            (P + 0.5*b_sqr)*delta(dir, 0) - b_con[dir]*b_cov[0];
//
//            flux[j][i][u1] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[1] + 
//                             (P + 0.5*b_sqr)*delta(dir, 1) - b_con[dir]*b_cov[1];
//
//            flux[j][i][u2] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[2] + 
//                             (P + 0.5*b_sqr)*delta(dir, 2) - b_con[dir]*b_cov[2];
//
//            flux[j][i][u3] = (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*u_cov[3] + 
//                             (P + 0.5*b_sqr)*delta(dir, 3) - b_con[dir]*b_cov[3];
//        }
//
//    return(0);
//
//}
//
//PetscScalar slope_lim(PetscScalar y1, PetscScalar y2, PetscScalar y3, 
//                      AppCtx *ctx)
//{
//    PetscScalar Dqm, Dqp, Dqc, s;
//
//    /* Monotonized Central or Woodward limiter */
//	if (ctx->lim == LIM_MC) {
//		Dqm = 2. * (y2 - y1);
//		Dqp = 2. * (y3 - y2);
//		Dqc = 0.5 * (y3 - y1);
//		s = Dqm * Dqp;
//		if (s <= 0.) {
//			return 0.;
//        }
//		else {
//			if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
//				return (Dqm);
//			else if (fabs(Dqp) < fabs(Dqc))
//				return (Dqp);
//			else
//				return (Dqc);
//		}
//	}
//	/* van leer slope limiter */
//	else if (ctx->lim == LIM_VANL) {
//		Dqm = (y2 - y1);
//		Dqp = (y3 - y2);
//		s = Dqm * Dqp;
//		if (s <= 0.)
//			return 0.;
//		else
//			return (2. * s / (Dqm + Dqp));
//	}
//
//	/* minmod slope limiter (crude but robust) */
//	else if (ctx->lim == LIM_MINM) {
//		Dqm = (y2 - y1);
//		Dqp = (y3 - y2);
//		s = Dqm * Dqp;
//		if (s <= 0.)
//			return 0.;
//		else if (fabs(Dqm) < fabs(Dqp))
//			return Dqm;
//		else
//			return Dqp;
//	}
//
//    SETERRQ(PETSC_COMM_SELF, 1, "Unknown dim slope limiter\n");
//    
//}
//
//PetscErrorCode Reconstruct2D(PetscScalar ***prim_lx, PetscScalar ***prim_rx,
//                             PetscScalar ***prim_ly, PetscScalar ***prim_ry,
//                             PetscScalar ***x, AppCtx *ctx)
//{
//    if (ctx->dim!=2) 
//        SETERRQ1(PETSC_COMM_SELF, 1, "2D Reconstruct called for %d dim\n", ctx->dim);
//
//    PetscInt i, j, var;
//    PetscScalar slope;
//
//    for (j=ctx->ystart-2; j<ctx->ystart + ctx->ysize+2; j++)
//        for (i=ctx->xstart-2; i<ctx->xstart + ctx->xsize+2; i++)
//            for (var=0; var<dof; var++) {
//                slope = slope_lim(x[j][i-1][var], x[j][i][var],
//                                  x[j][i+1][var], ctx);
//
//                prim_lx[j][i][var] = x[j][i][var] - 0.5*slope;
//                prim_rx[j][i][var] = x[j][i][var] + 0.5*slope;
//
//                slope = slope_lim(x[j-1][i][var], x[j][i][var],
//                                  x[j+1][i][var], ctx);
//                    
//                prim_ly[j][i][var] = x[j][i][var] - 0.5*slope;
//                prim_ry[j][i][var] = x[j][i][var] + 0.5*slope;
//
//            }
//
//    return(0);
//}
//
//PetscErrorCode ComputedU_dt2D(PetscScalar ***dU_dt, 
//                              PetscScalar ***x, PetscScalar ***dx_dt,
//                              AppCtx *ctx)
//{
//    PetscInt i, j, dir=0;
//    PetscScalar gamma, u_cov[4], u_con[4], b_cov[4], b_con[4], b_sqr, P;
//    PetscScalar dgamma_dt, du_cov_dt[4], du_con_dt[4];
//    PetscScalar db_cov_dt[4], db_con_dt[4], db_sqr_dt, dP_dt;
//
//    for (j=ctx->ystart; j<ctx->ystart + ctx->ysize; j++)
//        for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++) {
//            P = (5./3 - 1)*x[j][i][u];
//            dP_dt = (5./3 - 1)*dx_dt[j][i][u];
//
//            gamma = PetscSqrtScalar(1. + x[j][i][u1]*x[j][i][u1] + x[j][i][u2]*x[j][i][u2] +
//                                         x[j][i][u3]*x[j][i][u3]);
//            dgamma_dt = (x[j][i][u1]*dx_dt[j][i][u1] + x[j][i][u2]*dx_dt[j][i][u2] +
//                         x[j][i][u3]*dx_dt[j][i][u3])/gamma;
//
//            u_con[0] = gamma;
//            u_con[1] = x[j][i][u1];
//            u_con[2] = x[j][i][u2];
//            u_con[3] = x[j][i][u3];
//
//            du_con_dt[0] = dgamma_dt;
//            du_con_dt[1] = dx_dt[j][i][u1];
//            du_con_dt[2] = dx_dt[j][i][u2];
//            du_con_dt[3] = dx_dt[j][i][u3];
//
//            u_cov[0] = -gamma;
//            u_cov[1] = x[j][i][u1];
//            u_cov[2] = x[j][i][u2];
//            u_cov[3] = x[j][i][u3];
//
//            du_cov_dt[0] = -dgamma_dt;
//            du_cov_dt[1] = dx_dt[j][i][u1];
//            du_cov_dt[2] = dx_dt[j][i][u2];
//            du_cov_dt[3] = dx_dt[j][i][u3];
//
//            b_con[0] = x[j][i][B1]*u_con[1] + x[j][i][B2]*u_con[2] + x[j][i][B3]*u_con[3];
//            b_con[1] = (x[j][i][B1] + b_con[0]*u_con[1])/u_con[0];
//            b_con[2] = (x[j][i][B2] + b_con[0]*u_con[2])/u_con[0];
//            b_con[3] = (x[j][i][B3] + b_con[0]*u_con[3])/u_con[0];
//
//            db_con_dt[0] = dx_dt[j][i][B1]*u_con[1] + dx_dt[j][i][B2]*u_con[2] + 
//                           dx_dt[j][i][B3]*u_con[3] + x[j][i][B1]*du_con_dt[1] + 
//                           x[j][i][B2]*du_con_dt[2] + x[j][i][B3]*du_con_dt[3];
//
//            db_con_dt[1] = -(x[j][i][B1] + b_con[0]*u_con[1])*du_con_dt[0]/
//                            (u_con[0]*u_con[0]) + 
//                            (dx_dt[j][i][B1] + b_con[0]*du_con_dt[1] + 
//                             db_con_dt[0]*u_con[1])/u_con[0];
//        
//            db_con_dt[2] = -(x[j][i][B2] + b_con[0]*u_con[2])*du_con_dt[0]/
//                            (u_con[0]*u_con[0]) + 
//                            (dx_dt[j][i][B2] + b_con[0]*du_con_dt[2] + 
//                             db_con_dt[0]*u_con[2])/u_con[0];
//
//            db_con_dt[3] = -(x[j][i][B3] + b_con[0]*u_con[3])*du_con_dt[0]/
//                            (u_con[0]*u_con[0]) + 
//                            (dx_dt[j][i][B3] + b_con[0]*du_con_dt[3] + 
//                             db_con_dt[0]*u_con[3])/u_con[0];
//
//            b_cov[0] = -b_con[0];
//            b_cov[1] = b_con[1];
//            b_cov[2] = b_con[2];
//            b_cov[3] = b_con[3];
//
//            db_cov_dt[0] = -db_con_dt[0];
//            db_cov_dt[1] = db_con_dt[1];
//            db_cov_dt[2] = db_con_dt[2];
//            db_cov_dt[3] = db_con_dt[3];
//
//            b_sqr = -b_con[0]*b_con[0] + b_con[1]*b_con[1] + b_con[2]*b_con[2] +
//                     b_con[3]*b_con[3];
//
//            db_sqr_dt = 2.*(-b_con[0]*db_con_dt[0] + b_con[1]*db_con_dt[1] + 
//                            b_con[2]*db_con_dt[2] + b_con[3]*db_con_dt[3]);
//
//            dU_dt[j][i][rho] = dx_dt[j][i][rho]*u_con[dir] + x[j][i][rho]*du_con_dt[dir];
//
//            dU_dt[j][i][B1] = dx_dt[j][i][B1];
//
//            dU_dt[j][i][B2] = dx_dt[j][i][B2];
//
//            dU_dt[j][i][B3] = dx_dt[j][i][B3];
//
//
//            dU_dt[j][i][u] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
//                           db_sqr_dt)*u_con[dir]*u_cov[0] + 
//                          (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[0]
//                          + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[0]
//                          + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 0)
//                          - db_con_dt[dir]*b_cov[0] - b_con[dir]*db_cov_dt[0];
//
//            dU_dt[j][i][u1] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
//                            db_sqr_dt)*u_con[dir]*u_cov[1] + 
//                           (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[1]
//                           + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[1]
//                           + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 1)
//                           - db_con_dt[dir]*b_cov[1] - b_con[dir]*db_cov_dt[1];
//        
//            dU_dt[j][i][u2] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
//                            db_sqr_dt)*u_con[dir]*u_cov[2] + 
//                           (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[2]
//                           + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[2]
//                           + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 2)
//                           - db_con_dt[dir]*b_cov[2] - b_con[dir]*db_cov_dt[2];
//
//            dU_dt[j][i][u3] = (dP_dt + dx_dt[j][i][rho] + dx_dt[j][i][u] +
//                            db_sqr_dt)*u_con[dir]*u_cov[3] + 
//                           (P + x[j][i][rho] + x[j][i][u] + b_sqr)*du_con_dt[dir]*u_cov[3]
//                           + (P + x[j][i][rho] + x[j][i][u] + b_sqr)*u_con[dir]*du_cov_dt[3]
//                           + (dP_dt + 0.5*db_sqr_dt)*delta(dir, 3)
//                           - db_con_dt[dir]*b_cov[3] - b_con[dir]*db_cov_dt[3];
//
//    }
//    return(0);
//}
//
//
//PetscErrorCode SetBoundaryConditions2D(PetscScalar ***x,
//                                       PetscScalar ***dx_dt,
//                                       AppCtx *ctx)
//{
//    if (ctx->dim!=1) 
//        SETERRQ1(PETSC_COMM_SELF, 1, "2D SetBoundaryConditions called for %d dim\n", ctx->dim);
//
//    if (ctx->bdry_x == BDRY_PERIODIC)
//        continue;
//    if (ctx->bdry_x == BDRY_OUTFLOW) {
//        PetscInt i, j, var;
//
//        for (j=ctx->ystart; j<ctx->ystart+ctx->ysize; j++) {
//
//            // Inner r
//            for (i=ctx->xstart-ctx->nghost; i<ctx->xstart; i++)
//                for (var=0; var < dof; var++)
//                    x[j][i][var] = x[ctx->xstart][var];
//
//            // Outer r
//            for (i=ctx->xstart+ctx->xsize; 
//                 i<ctx->xstart+ctx->xsize+ctx->nghost; i++)
//                for (var=0; var < dof; var++)
//                    x[i][var] = x[ctx->xstart+ctx->xsize][var];
//
//        }
//
//    }
//    
//}
//
//PetscErrorCode FluxCT2D(PetscScalar ***flux_x,
//                        PetscScalar ***flux_y, AppCtx *ctx)
//{
//    if (ctx->dim!=2) 
//        SETERRQ1(PETSC_COMM_SELF, 1, "2D FluxCT called for %d dim\n", ctx->dim);
//
//    PetscInt i, j, i_local, j_local;
//    PetscScalar emf[ctx->ysize+ctx->nghost][ctx->xsize+ctx->nghost];
//
//    for (j=ctx->ystart; j<ctx->ystart + ctx->ysize + 2; j++)
//        for (i=ctx->xstart; i<ctx->xstart + ctx->xsize + 2; i++) {
//            i_local = i - ctx->xstart;
//            j_local = j - ctx->ystart;
//
//            emf[j_local][i_local] = 0.25*(flux_x[j][i][B2]  + 
//                                          flux_x[j-1][i][B2] -
//                                          flux_y[j][i][B1] - 
//                                          flux_y[j][i-1][B1]);
//        }
//
//    for (j=ctx->ystart; j<ctx->ystart + ctx->ysize + 1; j++)
//        for (i=ctx->xstart; i<ctx->xstart + ctx->xsize + 1; i++) {
//            i_local = i - ctx->xstart;
//            j_local = j - ctx->ystart;
//
//            flux_x[j][i][B1] = 0.;
//            
//            flux_x[j][i][B2] = 0.5*(emf[j_local][i_local]  +
//                                    emf[j_local+1][i_local]);
//
//            flux_y[j][i][B1] = -0.5*(emf[j_local][i_local] + 
//                                     emf[j_local][i_local+1]);
//
//            flux_y[j][i][B2] = 0.0;
//        }
//
//    return(0);
//
//}
//
//PetscErrorCode RiemannSolver2D(PetscScalar ***flux_x, 
//                               PetscScalar ***flux_y,
//                               PetscScalar ***flux_lx, 
//                               PetscScalar ***flux_rx,
//                               PetscScalar ***flux_ly, 
//                               PetscScalar ***flux_ry,
//                               PetscScalar ***U_lx, 
//                               PetscScalar ***U_rx,
//                               PetscScalar ***U_ly, 
//                               PetscScalar ***U_ry, AppCtx *ctx)
//{
//    if (ctx->dim!=2) 
//        SETERRQ1(PETSC_COMM_SELF, 1, "2D RiemannSolver called for %d dim\n", ctx->dim);
//
//    PetscInt i, j, var;
//
//    for (j=ctx->ystart-2; j<ctx->ystart + ctx->ysize+2; j++)
//        for (i=ctx->xstart-2; i<ctx->xstart + ctx->xsize+2; i++)
//            for (var=0; var<dof; var++) {
//                flux_x[j][i][var] = 0.5*(flux_lx[j][i][var] + flux_rx[j][i-1][var] -
//                                        (U_lx[j][i][var] - U_rx[j][i-1][var]));
//                flux_y[j][i][var] = 0.5*(flux_ly[j][i][var] + flux_ry[j-1][i][var] -
//                                        (U_ly[j][i][var] - U_ry[j-1][i][var]));
//            }
//
//    return(0);
//
//}
//
//
//PetscErrorCode ResFunction(TS ts,
//                           PetscScalar t, 
//                           Vec Xvec, Vec dX_dt, Vec F,
//                           void *ptr)
//{
//    AppCtx *ctx = (AppCtx*)ptr;
//    
//    DM da;
//
//    Vec localX;
//
//    TSGetDM(ts, &da);
//
//    DMGetLocalVector(da, &localX);
//
//    DMDAGetCorners(da, 
//                   &ctx->xstart, &ctx->ystart, &ctx->zstart,
//                   &ctx->xsize, &ctx->ysize, &ctx->zsize);
//
//    DMGlobalToLocalBegin(da, Xvec, INSERT_VALUES, localX);
//    DMGlobalToLocalEnd(da, Xvec, INSERT_VALUES, localX);
//
//    switch (ctx->dim) {
//        case 1:
//        {
//
//            break;
//        }
//        case 2:
//        {
//            PetscInt i, j, var;
//            PetscScalar ***x, ***f, ***dx_dt;
//
//            DMDAVecGetArrayDOF(da, localX, &x);
//            DMDAVecGetArrayDOF(da, F, &f);
//            DMDAVecGetArrayDOF(da, dX_dt, &dx_dt);
//
//            SetBoundaryConditions2D(x, dx_dt, ctx);
//
//            ComputedU_dt2D(f, x, dx_dt, ctx);
//
//            Reconstruct2D(prim_lx, prim_rx,
//                          prim_ly, prim_ry,
//                          x, ctx);
//
//            ComputeFlux2D(flux_lx, prim_lx, 1, ctx);
//            ComputeFlux2D(flux_rx, prim_rx, 1, ctx);
//            ComputeFlux2D(flux_ly, prim_ly, 2, ctx);
//            ComputeFlux2D(flux_ry, prim_ry, 2, ctx);
//
//            ComputeFlux2D(U_lx, prim_lx, 0, ctx);
//            ComputeFlux2D(U_rx, prim_rx, 0, ctx);
//            ComputeFlux2D(U_ly, prim_ly, 0, ctx);
//            ComputeFlux2D(U_ry, prim_ry, 0, ctx);
//
//            RiemannSolver2D(flux_x, flux_y,
//                            flux_lx, flux_rx,
//                            flux_ly, flux_ry,
//                            U_lx, U_rx, 
//                            U_ly, U_ry, ctx);
//
//            FluxCT2D(flux_x, flux_y, ctx);
//
//            for (j=ctx->ystart; j<ctx->ystart + ctx->ysize; j++)
//                for (i=ctx->xstart; i<ctx->xstart + ctx->xsize; i++)
//                    for (var=0; var<dof; var++) {
//                        f[j][i][var] =  f[j][i][var] + 
//                                        (flux_x[j][i+1][var] - 
//                                         flux_x[j][i][var])/ctx->dx +
//                                        (flux_y[j+1][i][var] -
//                                         flux_y[j][i][var])/ctx->dy;
//                    }
//
//
//            DMDAVecRestoreArrayDOF(da, localX, &x);
//            DMDAVecRestoreArrayDOF(da, F, &f);
//            DMDAVecRestoreArrayDOF(da, dX_dt, &dx_dt);
//
//            break;
//        }
//        case 3:
//        {
//
//            break;
//        }
//
//        default: SETERRQ(PETSC_COMM_SELF, 1, "Unknown dim in RHSFunction\n");
//
//    }
//
//    DMRestoreLocalVector(da, &localX);
//    return(0);
//}
