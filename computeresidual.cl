#include "constants.h"
#include "reconstruct.cl"
#include "phys.cl"
#include "riemannsolver.cl"

__kernel void ComputeResidual(__global const REAL* restrict prim, 
                              __global const REAL* restrict dprim_dt,
                              __global REAL* restrict F)
{
    int i = get_global_id(0);
    int j = get_global_id(1);
    int iTile = get_local_id(0);
    int jTile = get_local_id(1);
    
    // Finite Volume variables
    REAL primEdge[DOF];
    REAL fluxL[DOF], uL[DOF], fluxR[DOF], uR[DOF];
    REAL fluxX1L[DOF], fluxX1R[DOF], fluxX2L[DOF], fluxX2R[DOF];
    REAL cminL, cminR, cmaxL, cmaxR;

    // Geometry variables
    REAL X1, X2;
    REAL gcon[NDIM][NDIM], gcov[NDIM][NDIM], gdet;
    REAL alpha, gamma, g;

    // Physics variables
    REAL dU_dt[DOF], dvars_dt[DOF];
    REAL ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
    REAL mhd[NDIM][NDIM], bsqr;

    // Time derivative of physics variables
    REAL ducon_dt[NDIM], ducov_dt[NDIM], dbcon_dt[NDIM], dbcov_dt[NDIM];
    REAL dgamma_dt;

    __local REAL primTile[(TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*DOF];

    for (int var=0; var<DOF; var++) {
        primTile[INDEX_LOCAL(iTile,jTile,var)] = prim[INDEX_GLOBAL(i,j,var)];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

//    primEdge[RHO] = exp(primTile[INDEX_LOCAL(iTile,jTile,RHO)]);
//    dvars_dt[RHO] = primEdge[RHO]*dprim_dt[INDEX_GLOBAL(i,j,RHO)];
//
//    primEdge[UU] = exp(primTile[INDEX_LOCAL(iTile,jTile,UU)]);
//    dvars_dt[UU] = primEdge[UU]*dprim_dt[INDEX_GLOBAL(i,j,UU)];
//
//    if (i==10 && j==30)
//    printf("i = %d, j = %d, primEdge[RHO] = %.15f\n", i, j, primEdge[RHO]);

    for (int var=0; var<DOF; var++) {
        primEdge[var] = primTile[INDEX_LOCAL(iTile,jTile,var)];
        dvars_dt[var] = dprim_dt[INDEX_GLOBAL(i,j,var)];
    }
    

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    dgammaCalc_dt(&dgamma_dt, gamma, primEdge, dvars_dt, gcov);
    duconCalc_dt(ducon_dt, dgamma_dt, alpha, dvars_dt, gcon);
    covFromCon(ducov_dt, ducon_dt, gcov);

    dbconCalc_dt(dbcon_dt, ucon, ducon_dt, ucov, ducov_dt, bcon, 
                 primEdge, dvars_dt);
    covFromCon(dbcov_dt, dbcon_dt, gcov);

    ComputedU_dt(dU_dt,
                 ucon, ducon_dt, ucov, ducov_dt,
                 bcon, dbcon_dt, bcov, dbcov_dt,
                 gcon, gcov, primEdge, dvars_dt,
                 gamma, dgamma_dt, alpha, g);

    addSources(dU_dt, 
               ucon, ucov, bcon, bcov, 
               gcon, gcov, mhd, primEdge, g,
               X1, X2);

    for (int var=0; var<DOF; var++) {

        if (iTile==0) {
            if (i>=NG) {
                for (int iNg=-NG; iNg<0; iNg++) {
                    primTile[INDEX_LOCAL(iNg,jTile,var)] = 
                        prim[INDEX_GLOBAL(i+iNg,j,var)];
                }
            } else {
                for (int iNg=-NG; iNg<0; iNg++) {
                    primTile[INDEX_LOCAL(iNg,jTile,var)] = 
                    primTile[INDEX_LOCAL(0,jTile,var)];
                }
//                for (int iNg=-NG; iNg<0; iNg++) {
//                    primTile[INDEX_LOCAL(iNg,jTile,var)] =
//                        prim[INDEX_GLOBAL(N1+iNg,j,var)];
//                }
            }
        }
    
        if (iTile==TILE_SIZE_X1-1) {
            if (i<=N1-NG) {
                for (int iNg=0; iNg<NG; iNg++) {
                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] = 
                        prim[INDEX_GLOBAL(i+iNg+1,j,var)];
                }
            } else {
                for (int iNg=0; iNg<NG; iNg++) {
                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] =
                    primTile[INDEX_LOCAL(TILE_SIZE_X1-1,jTile,var)];
                }
//                for (int iNg=0; iNg<NG; iNg++) {
//                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] =
//                    prim[INDEX_GLOBAL(iNg,j,var)];
//                }
            }
        }
       
        if (jTile==0) {
            if (j>=NG) {
                for (int jNg=-NG; jNg<0; jNg++) {
                    primTile[INDEX_LOCAL(iTile,jNg,var)] = 
                        prim[INDEX_GLOBAL(i,j+jNg,var)];
                }
            } else {
                for (int jNg=-NG; jNg<0; jNg++) {
                    primTile[INDEX_LOCAL(iTile,jNg,var)] = 
                    primTile[INDEX_LOCAL(iTile,-jNg-1,var)];
                }
//                for (int jNg=-NG; jNg<0; jNg++) {
//                    primTile[INDEX_LOCAL(iTile,jNg,var)] = 
//                    prim[INDEX_GLOBAL(i,N2+jNg,var)];
//                }
            }
        }
    
        if (jTile==TILE_SIZE_X2-1) {
            if (j<=N2-NG) {
                for (int jNg=0; jNg<NG; jNg++) {
                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] = 
                        prim[INDEX_GLOBAL(i,j+jNg+1,var)];
                }
            } else {
                for (int jNg=0; jNg<NG; jNg++) {
                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] =
                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2-1-jNg,var)];
                }
//                for (int jNg=0; jNg<NG; jNg++) {
//                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] =
//                        prim[INDEX_GLOBAL(i,jNg,var)];
//                }
            }
        }
    
    }

    if (jTile==0 && j<NG) {
        for (int jNg=-NG; jNg<0; jNg++) {
            primTile[INDEX_LOCAL(iTile,jNg,U2)] = 
           -primTile[INDEX_LOCAL(iTile,jNg,U2)];
            primTile[INDEX_LOCAL(iTile,jNg,B2)] = 
           -primTile[INDEX_LOCAL(iTile,jNg,B2)];
        }
    }

    if (jTile==TILE_SIZE_X2-1 && j>N2-NG) {
        for (int jNg=0; jNg<NG; jNg++) {
            primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,U2)] = 
           -primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,U2)];
            primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,B2)] = 
           -primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,B2)];
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    
/*
     (i,j) 
     _____
    |     |
    |  o  | X1 = i_TO_X1_CENTER(i)
    |_____| X2 = j_TO_X2_CENTER(j)
     _____
    |     |
    |o    | X1 = i_TO_X1_FACE(i)
    |_____| X2 = j_TO_X2_CENTER(j)

     _____
    |     |
    |     | X1 = i_TO_X1_CENTER(i)
    |__o__| X2 = j_TO_X2_FACE(j)

     _____
    |     |
    |     | X1 = i_TO_X1_FACE(i)
    |o____| X2 = j_TO_X2_FACE(j)

*/

/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |  _  |  _  |  _  |     |     |
    |  o  |  o->|  o  |  o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile-1, jTile, primEdge, RIGHT);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_FACE(i); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 1); 

/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |     |  _  |  _  |  _  |     |
    |  o  |  o  |<-o  |  o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile, jTile, primEdge, LEFT);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_FACE(i); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 1); 

    RiemannSolver(fluxX1L, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);


/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |     |  _  |  _  |  _  |     |
    |  o  |  o  |  o->|  o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile, jTile, primEdge, RIGHT);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_FACE(i+1); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 1); 

/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |     |     |  _  |  _  |  _  |
    |  o  |  o  |  o  |<-o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile+1, jTile, primEdge, LEFT);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_FACE(i+1); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 1); 

    RiemannSolver(fluxX1R, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);


/*
                 _____ 
                |     |
                |  o  |  j+2
                |_____| 
                |     |
                |  o  |  j+1 
                |_____| 
                |     |
                | |o| |   j
                |_____|
                |  ^  |
                | |o| |  j-1
                |_____| 
                |     |
                | |o| |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile-1, primEdge, UP);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 2); 

/*
                 _____ 
                |     |
                |  o  |  j+2
                |_____| 
                |     |
                | |o| |  j+1 
                |_____| 
                |     |
                | |o| |   j
                |__v__|
                |     |
                | |o| |  j-1
                |_____| 
                |     |
                |  o  |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile, primEdge, DOWN);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 2); 

    RiemannSolver(fluxX2L, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);


/*
                 _____ 
                |     |
                |  o  |  j+2
                |_____| 
                |     |
                | |o| |  j+1 
                |_____| 
                |  ^  |
                | |o| |   j
                |_____|
                |     |
                | |o| |  j-1
                |_____| 
                |     |
                |  o  |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile, primEdge, UP);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j+1);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 2); 

/*
                 _____ 
                |     |
                | |o| |  j+2
                |_____| 
                |     |
                | |o| |  j+1 
                |__V__| 
                |     |
                | |o| |   j
                |_____|
                |     |
                |  o  |  j-1
                |_____| 
                |     |
                |  o  |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile+1, primEdge, DOWN);

//    primEdge[RHO] = exp(primEdge[RHO]);
//    primEdge[UU] = exp(primEdge[UU]);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j+1);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 2); 

    RiemannSolver(fluxX2R, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);

    for (int var=0; var<DOF; var++) {
        F[INDEX_GLOBAL(i,j,var)] = (dU_dt[var] +
                                  (fluxX1R[var] - fluxX1L[var])/DX1 +
                                  (fluxX2R[var] - fluxX2L[var])/DX2);
    }

}
























//__kernel void ComputeResidual(__global const REAL* restrict prim, 
//                              __global const REAL* restrict dprim_dt,
//                              __global REAL* restrict F)
//{
//
//    int i = get_global_id(0);
//    int j = get_global_id(1);
//    int iTile = get_local_id(0);
//    int jTile = get_local_id(1);
//
//    __local REAL primTile[(TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*DOF];
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
//    REAL mhd[NDIM][NDIM], vars[DOF], dvars_dt[DOF];
//
//    for (int var=0; var<DOF; var++) {
//        primTile[INDEX_LOCAL(iTile,jTile,var)] = prim[INDEX_GLOBAL(i,j,var)];
//    }
//
//    barrier(CLK_LOCAL_MEM_FENCE);
//
//    for (int var=0; var<DOF; var++) {
//        vars[var] = primTile[INDEX_LOCAL(iTile,jTile,var)];
//        dvars_dt[var] = dprim_dt[INDEX_GLOBAL(i,j,var)];
//    }
//
//    for (int var=0; var<DOF; var++) {
//        vars[var] = prim[INDEX_GLOBAL(i,j,var)];
//    }
//
//    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j);
////    ComputedU_dt(vars, dvars_dt, dU_dt, X1, X2);
//    for (int n=0; n<DOF; n++)
//        dU_dt[n] = 0.;
//
////    addSources(dU_dt,
////               ucon, ucov, bcon, bcov,
////               gcon, gcov, mhd, vars, 
////               X1, X2);
//
//    for (int var=0; var<DOF; var++) {
//
//        if (iTile==0) {
//            if (i>=NG) {
//                for (int iNg=-NG; iNg<0; iNg++) {
//                    primTile[INDEX_LOCAL(iNg,jTile,var)] = 
//                        prim[INDEX_GLOBAL(i+iNg,j,var)];
//                }
//            } else {
//                for (int iNg=-NG; iNg<0; iNg++) {
//                    primTile[INDEX_LOCAL(iNg,jTile,var)] = 
//                    primTile[INDEX_LOCAL(0,jTile,var)];
//                }
////                for (int iNg=-NG; iNg<0; iNg++) {
////                    primTile[INDEX_LOCAL(iNg,jTile,var)] =
////                        prim[INDEX_GLOBAL(N1+iNg,j,var)];
////                }
//            }
//        }
//    
//        if (iTile==TILE_SIZE_X1-1) {
//            if (i<=N1-NG) {
//                for (int iNg=0; iNg<NG; iNg++) {
//                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] = 
//                        prim[INDEX_GLOBAL(i+iNg+1,j,var)];
//                }
//            } else {
//                for (int iNg=0; iNg<NG; iNg++) {
//                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] =
//                    primTile[INDEX_LOCAL(TILE_SIZE_X1-1,jTile,var)];
//                }
////                for (int iNg=0; iNg<NG; iNg++) {
////                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] =
////                    prim[INDEX_GLOBAL(iNg,j,var)];
////                }
//            }
//        }
//       
//        if (jTile==0) {
//            if (j>=NG) {
//                for (int jNg=-NG; jNg<0; jNg++) {
//                    primTile[INDEX_LOCAL(iTile,jNg,var)] = 
//                        prim[INDEX_GLOBAL(i,j+jNg,var)];
//                }
//            } else {
//                for (int jNg=-NG; jNg<0; jNg++) {
//                    primTile[INDEX_LOCAL(iTile,jNg,var)] = 
//                    primTile[INDEX_LOCAL(iTile,-jNg-1,var)];
//                }
////                for (int jNg=-NG; jNg<0; jNg++) {
////                    primTile[INDEX_LOCAL(iTile,jNg,var)] = 
////                    prim[INDEX_GLOBAL(i,N2+jNg,var)];
////                }
//            }
//        }
//    
//        if (jTile==TILE_SIZE_X2-1) {
//            if (j<=N2-NG) {
//                for (int jNg=0; jNg<NG; jNg++) {
//                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] = 
//                        prim[INDEX_GLOBAL(i,j+jNg+1,var)];
//                }
//            } else {
//                for (int jNg=0; jNg<NG; jNg++) {
//                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] =
//                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2-1-jNg,var)];
//                }
////                for (int jNg=0; jNg<NG; jNg++) {
////                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] =
////                        prim[INDEX_GLOBAL(i,jNg,var)];
////                }
//            }
//        }
//    
//    }
//
//    if (jTile==0 && j<NG) {
//        for (int jNg=-NG; jNg<0; jNg++) {
//            primTile[INDEX_LOCAL(iTile,jNg,U2)] = 
//           -primTile[INDEX_LOCAL(iTile,jNg,U2)];
//            primTile[INDEX_LOCAL(iTile,jNg,B2)] = 
//           -primTile[INDEX_LOCAL(iTile,jNg,B2)];
//        }
//    }
//
//    if (jTile==TILE_SIZE_X2-1 && j>N2-NG) {
//        for (int jNg=0; jNg<NG; jNg++) {
//            primTile[INDEX_LOCAL(iTile,jNg,U2)] = 
//           -primTile[INDEX_LOCAL(iTile,jNg,U2)];
//            primTile[INDEX_LOCAL(iTile,jNg,B2)] = 
//           -primTile[INDEX_LOCAL(iTile,jNg,B2)];
//        }
//    }
//
//    barrier(CLK_LOCAL_MEM_FENCE);
//
//    // Compute fluxes along X1
//    X1 = i_TO_X1_FACE(i); X2 = j_TO_X2_CENTER(j);
//    ReconstructX1(primTile, iTile-1, jTile, primEdge, RIGHT);
//    ComputeFluxAndU(fluxR, uR, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 1);
//
//    ReconstructX1(primTile, iTile, jTile, primEdge, LEFT);
//    ComputeFluxAndU(fluxL, uL, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 1);
//
//    RiemannSolver(fluxR, fluxL, uR, uL, fluxX1L);
//
//
//    X1 = i_TO_X1_FACE(i+1); X2 = j_TO_X2_CENTER(j);
//    ReconstructX1(primTile, iTile, jTile, primEdge, RIGHT);
//    ComputeFluxAndU(fluxR, uR, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 1);
//
//    ReconstructX1(primTile, iTile+1, jTile, primEdge, LEFT);
//    ComputeFluxAndU(fluxL, uL, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 1);
//
//    RiemannSolver(fluxR, fluxL, uR, uL, fluxX1R);
//
//
//    // Compute fluxes along X2
//    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j);
//    ReconstructX2(primTile, iTile, jTile-1, primEdge, RIGHT);
//    ComputeFluxAndU(fluxR, uR, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 2);
//
//    ReconstructX2(primTile, iTile, jTile, primEdge, LEFT);
//    ComputeFluxAndU(fluxL, uL, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 2);
//
//    RiemannSolver(fluxR, fluxL, uR, uL, fluxX2L);
//
//
//    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j+1);
//    ReconstructX2(primTile, iTile, jTile, primEdge, RIGHT);
//    ComputeFluxAndU(fluxR, uR, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 2);
//
//    ReconstructX2(primTile, iTile, jTile+1, primEdge, LEFT);
//    ComputeFluxAndU(fluxL, uL, 
//                    ucon, ucov, bcon, bcov, 
//                    gcon, gcov, mhd, primEdge, 
//                    X1, X2, 2);
//
//    RiemannSolver(fluxR, fluxL, uR, uL, fluxX2R);
//
//
//    for (int var=0; var<DOF; var++) {
//        F[INDEX_GLOBAL(i,j,var)] = dU_dt[var] +
//                                   (fluxX1R[var] - fluxX1L[var])/DX1 + 
//                                   (fluxX2R[var] - fluxX2L[var])/DX2;
//
////        if (isnan(F[INDEX_GLOBAL(i,j,var)]))
////        printf("i = %d, j = %d, iTile = %d, jTile = %d, var = %d, vars = %f, f = %f\n",
////                i, j, iTile, jTile, var, vars[var], F[INDEX_GLOBAL(i,j,var)]);
//    }
//}



