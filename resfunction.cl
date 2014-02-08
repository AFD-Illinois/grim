#include "constants.h"

#define INDEX(i,j,varNum) (i+1+(TILE_SIZE+2)*(((j)%(2*NG+1))+(2*NG+1)*varNum))

#define LOAD_VAR(i, j, iOffset, jOffset, array, varNum) \
    (VLOAD_SIZE(0, \
                (DOUBLE*)(&array[INDEX(i, j + jOffset, varNum)]) + iOffset) )

#define STORE_VAR(val, i, j, iOffset, jOffset, array, varNum) \
    (VSTORE_SIZE(val, 0, \
                (DOUBLE*)(&array[INDEX(i, j + jOffset, varNum)]) + iOffset) )

//void RiemannSolver(const DOUBLE_SIZE* restrict flux_lx, 
//                   const DOUBLE_SIZE* restrict flux_rx,
//                   const DOUBLE_SIZE* restrict flux_ly,
//                   const DOUBLE_SIZE* restrict flux_ry,
//                   const DOUBLE_SIZE* restrict U_lx,
//                   const DOUBLE_SIZE* restrict U_rx,
//                   const DOUBLE_SIZE* restrict U_ly, 
//                   const DOUBLE_SIZE* restrict U_ry,
//                   const int iStart, const int iEnd,
//                   const int jStart, const int jEnd,
//                   DOUBLE_SIZE* restrict flux_x,
//                   DOUBLE_SIZE* restrict flux_y)
//{
//    DOUBLE_SIZE tmp;
//    for (int var=0;  var<DOF; var++)
//        for (int j=jStart+1; j<jEnd-1; j++) {
//            /*i=-1*/
//            tmp = LOAD_VAR(iStart, j, 1, 0, flux_lx, var) +
//                  LOAD_VAR(iStart, j, 0, 0, flux_rx, var) +
//                  LOAD_VAR(iStart, j, 1, 0, U_lx, var) -
//                  LOAD_VAR(iStart, j, 0, 0, U_rx, var);
//
//            STORE_VAR(tmp, iStart, j, 1, 0, flux_x, var);
//
//            tmp = LOAD_VAR(iStart, j, 0, 0, flux_ly, var) +
//                  LOAD_VAR(iStart, j, 0, -1, flux_ry, var) +
//                  LOAD_VAR(iStart, j, 0, 0, U_ly, var) -
//                  LOAD_VAR(iStart, j, 0, -1, U_ry, var);
//
//            STORE_VAR(tmp, iStart, j, 0, 0, flux_y, var);
//
//            for (int i=iStart+1; i<iEnd; i++) {
//                tmp = LOAD_VAR(i, j, 0, 0, flux_lx, var) +
//                      LOAD_VAR(i, j, -1, 0, flux_rx, var) +
//                      LOAD_VAR(i, j, 0, 0, U_lx, var) -
//                      LOAD_VAR(i, j, -1, 0, U_rx, var);
//
//                STORE_VAR(tmp, i, j, 0, 0, flux_x, var);
//
//                tmp = LOAD_VAR(i, j, 0, 0, flux_ly, var) +
//                      LOAD_VAR(i, j, 0, -1, flux_ry, var) +
//                      LOAD_VAR(i, j, 0, 0, U_ly, var) -
//                      LOAD_VAR(i, j, 0, -1, U_ry, var);
//
//                STORE_VAR(tmp, i, j, 0, 0, flux_y, var);
//            }
//        }
//
//}
//
//void FluxCT(const int iStart, const int iEnd,
//            const int jStart, const int jEnd,
//            DOUBLE_SIZE* restrict flux_x,
//            DOUBLE_SIZE* restrict flux_y)
//{
//    int iLocal, jLocal, iLocalSize, jLocalSize;
//    iLocalSize = iEnd - iStart;
//    jLocalSize = jEnd - jStart;
//    DOUBLE_SIZE emf[jLocalSize][iLocalSize], tmp;
//
////    for (int j=jStart+1; j<jEnd; j++)
////        for (int i=iStart+1; i<iEnd; i++) {
////            iLocal = i - iStart - 1;
////            jLocal = j - jStart - 1;
////
////            emf[jLocal][iLocal] =
////                        0.25*(LOAD_VAR(i, j, iOffset, 0, flux_x, B2) +
////                              LOAD_VAR(i, j, iOffset, -1, flux_x, B2) +
////                              LOAD_VAR(i, j, iOffset, 0, flux_y, B1) +
////                              LOAD_VAR(i, j, iOffset-1, 0, flux_y, B1));
////
////        }
////
////    for (int j=jStart+1; j<jEnd-1; j++)
////        for (int i=iStart+1; i<iEnd-1; i++) {
////            iLocal = i - iStart - 1;
////            jLocal = j - jStart - 1;
////
////            STORE_VAR(0., i, j, iOffset, 0, flux_x, B1);
////
////            tmp = 0.5*(emf[jLocal][iLocal] + emf[jLocal+1][iLocal]);
////
////            STORE_VAR(tmp, i, j, iOffset, 0, flux_x, B2);
////
////            tmp = VLOAD_SIZE(0, 
////                            (DOUBLE*)(&emf[iLocal + jLocal*iLocalSize])+1) -
////                  0.5*emf[jLocal][iLocal];
////
////            STORE_VAR(tmp, i, j, iOffset, 0, flux_y, B1);
////
////            STORE_VAR(0., i, j, iOffset, 0, flux_y, B2);
////        }
//
//}
//
//DOUBLE_SIZE SlopeLim(DOUBLE_SIZE y1, DOUBLE_SIZE y2, DOUBLE_SIZE y3)
//{
//    return 0.5*(y3 - y1);
//}
//
//void Reconstruct(const DOUBLE_SIZE* restrict prim, 
//                 const int iStart, const int iEnd,
//                 const int jStart, const int jEnd,
//                 DOUBLE_SIZE* restrict prim_lx,
//                 DOUBLE_SIZE* restrict prim_rx,
//                 DOUBLE_SIZE* restrict prim_ly, 
//                 DOUBLE_SIZE* restrict prim_ry)
//{
//    int i;
//    DOUBLE_SIZE mid, slope;
//
//    for (int var=0; var<DOF; var++)
//        for (int j=jStart+1; j<jEnd-1; j++) {
//        
//            /*i = -1*/
//            mid = LOAD_VAR(iStart, j, 1, 0, prim, var);
//
//            slope = SlopeLim(LOAD_VAR(iStart, j, 0, 0, prim, var),
//                             mid,
//                             LOAD_VAR(iStart, j, 2, 0, prim, var));
//
//            STORE_VAR(mid - 0.5*slope, iStart, j, 1, 0, prim_lx, var);
//            STORE_VAR(mid + 0.5*slope, iStart, j, 1, 0, prim_rx, var);
//
//            slope = SlopeLim(LOAD_VAR(iStart, j, 1, -1, prim, var),
//                             mid,
//                             LOAD_VAR(iStart, j, 1, 1, prim, var));
//
//            STORE_VAR(mid - 0.5*slope, iStart, j, 1, 0, prim_ly, var);
//            STORE_VAR(mid + 0.5*slope, iStart, j, 1, 0, prim_ry, var);
//
//            for (int i=iStart+1; i< iEnd-1; i++) {
//                mid = LOAD_VAR(i, j, 0, 0, prim, var);
//
//                slope = SlopeLim(LOAD_VAR(i, j, -1, 0, prim, var),
//                                 mid,
//                                 LOAD_VAR(i, j, 1, 0, prim, var));
//
//                STORE_VAR(mid - 0.5*slope, i, j, 0, 0, prim_lx, var);
//                STORE_VAR(mid + 0.5*slope, i, j, 0, 0, prim_rx, var);
//
//                slope = SlopeLim(LOAD_VAR(i, j, 0, -1, prim, var),
//                                 mid,
//                                 LOAD_VAR(i, j, 0, 1, prim, var));
//
//                STORE_VAR(mid - 0.5*slope, i, j, 0, 0, prim_ly, var);
//                STORE_VAR(mid + 0.5*slope, i, j, 0, 0, prim_ry, var);
//            }
//
//            /*i=TILE_SIZE*/
//            mid = LOAD_VAR(iEnd-1, j, -1, 0, prim, var);
//
//            slope = SlopeLim(LOAD_VAR(iEnd-1, j, -2, 0, prim, var),
//                             mid,
//                             LOAD_VAR(iEnd-1, j, 0, 0, prim, var));
//
//            STORE_VAR(mid - 0.5*slope, iEnd-1, j, -1, 0, prim_lx, var);
//            STORE_VAR(mid + 0.5*slope, iEnd-1, j, -1, 0, prim_rx, var);
//
//            slope = SlopeLim(LOAD_VAR(iEnd-1, j, -1, -1, prim, var),
//                             mid,
//                             LOAD_VAR(iEnd-1, j, -1, 1, prim, var));
//
//            STORE_VAR(mid - 0.5*slope, iEnd-1, j, -1, 0, prim_ly, var);
//            STORE_VAR(mid + 0.5*slope, iEnd-1, j, -1, 0, prim_ry, var);
//        }
//
//}
//
void ComputeFluxAndU(const DOUBLE_SIZE* restrict prim,
                     const int dir, 
                     const int iStart, const int iEnd,
                     const int jStart, const int jEnd,
                     DOUBLE_SIZE* restrict flux,
                     DOUBLE_SIZE* restrict U)
{

#define DELTA(p, q) (p==q ? 1. : 0.)
#define ALPHA (1.)
#define GCOV(mu, nu) (mu==0 ? -DELTA(mu, nu) : DELTA(mu, nu))

    DOUBLE_SIZE gamma, P, tmp1, tmp2, bsqr;
    DOUBLE_SIZE ucon[NDIM], ucov[NDIM];
    DOUBLE_SIZE bcon[NDIM], bcov[NDIM];
    DOUBLE_SIZE var[DOF];

    for (int j=jStart; j<jEnd; j++) {
        for (int i=iStart; i<iEnd; i++) {
        
            var[RHO] = prim[INDEX(i, j, RHO)];
            var[UU] = prim[INDEX(i, j, UU)];
            var[U1] = prim[INDEX(i, j, U1)];
            var[U2] = prim[INDEX(i, j, U2)];
            var[U3] = prim[INDEX(i, j, U3)];
            var[B1] = prim[INDEX(i, j, B1)];
            var[B2] = prim[INDEX(i, j, B2)];
            var[B3] = prim[INDEX(i, j, B3)];

            gamma = sqrt(1 + GCOV(1, 1)*var[U1]*var[U1] + 
                             GCOV(2, 2)*var[U2]*var[U2] +
                             GCOV(3, 3)*var[U3]*var[U3] +
                          2*(GCOV(1, 2)*var[U1]*var[U2] +
                             GCOV(1, 3)*var[U1]*var[U3] +
                             GCOV(2, 3)*var[U2]*var[U3]));

            ucon[0] = gamma/ALPHA;
            ucon[1] = var[U1] - gamma*GCOV(0, 1)*ALPHA;
            ucon[2] = var[U2] - gamma*GCOV(0, 2)*ALPHA;
            ucon[3] = var[U3] - gamma*GCOV(0, 3)*ALPHA;

            ucov[0] = GCOV(0, 0)*ucon[0] +
                      GCOV(0, 1)*ucon[1] + 
                      GCOV(0, 2)*ucon[2] + 
                      GCOV(0, 3)*ucon[3];

            ucov[1] = GCOV(1, 0)*ucon[0] + 
                      GCOV(1, 1)*ucon[1] + 
                      GCOV(1, 2)*ucon[2] + 
                      GCOV(1, 3)*ucon[3];

            ucov[2] = GCOV(2, 0)*ucon[0] + 
                      GCOV(2, 1)*ucon[1] + 
                      GCOV(2, 2)*ucon[2] + 
                      GCOV(2, 3)*ucon[3];

            ucov[3] = GCOV(3, 0)*ucon[0] + 
                      GCOV(3, 1)*ucon[1] + 
                      GCOV(3, 2)*ucon[2] + 
                      GCOV(3, 3)*ucon[3];

            bcon[0] = var[B1]*ucon[1] + var[B2]*ucon[2] + var[B3]*ucon[3];
            bcon[1] = (var[B1] + bcon[0]*ucon[1])/ucon[0];
            bcon[2] = (var[B2] + bcon[0]*ucon[2])/ucon[0];
            bcon[3] = (var[B3] + bcon[0]*ucon[3])/ucon[0];

            bcov[0] = GCOV(0, 0)*bcon[0] +
                      GCOV(0, 1)*bcon[1] + 
                      GCOV(0, 2)*bcon[2] + 
                      GCOV(0, 3)*bcon[3];

            bcov[1] = GCOV(1, 0)*bcon[0] + 
                      GCOV(1, 1)*bcon[1] + 
                      GCOV(1, 2)*bcon[2] + 
                      GCOV(1, 3)*bcon[3];

            bcov[2] = GCOV(2, 0)*bcon[0] + 
                      GCOV(2, 1)*bcon[1] + 
                      GCOV(2, 2)*bcon[2] + 
                      GCOV(2, 3)*bcon[3];

            bcov[3] = GCOV(3, 0)*bcon[0] + 
                      GCOV(3, 1)*bcon[1] + 
                      GCOV(3, 2)*bcon[2] + 
                      GCOV(3, 3)*bcon[3];

            bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] + 
                   bcon[2]*bcov[2] + bcon[3]*bcov[3];

            P = (GAMMA - 1.)*var[UU];
            tmp1 = P + var[RHO] + var[UU] + bsqr;
            tmp2 = P + 0.5*bsqr;

            flux[INDEX(i, j, RHO)]  = var[RHO]*ucon[dir];
            flux[INDEX(i, j, UU)] = tmp1*ucon[dir]*ucov[0] + tmp2*DELTA(dir, 0)
                                    - bcon[dir]*bcov[0];
            flux[INDEX(i, j, U1)] = tmp1*ucon[dir]*ucov[1] + tmp2*DELTA(dir, 1) 
                                    - bcon[dir]*bcov[1];
            flux[INDEX(i, j, U2)] = tmp1*ucon[dir]*ucov[2] + tmp2*DELTA(dir, 2) 
                                    - bcon[dir]*bcov[2];
            flux[INDEX(i, j, U3)] = tmp1*ucon[dir]*ucov[3] + tmp2*DELTA(dir, 3) 
                                    - bcon[dir]*bcov[3];

            flux[INDEX(i, j, B1)] = bcon[1]*ucon[dir] - bcon[dir]*ucon[1];
            flux[INDEX(i, j, B2)] = bcon[2]*ucon[dir] - bcon[dir]*ucon[2];
            flux[INDEX(i, j, B3)] = bcon[3]*ucon[dir] - bcon[dir]*ucon[3];

            U[INDEX(i, j, RHO)]  = var[RHO]*ucon[0];

            U[INDEX(i, j, B1)] = bcon[1]*ucon[0] - bcon[0]*ucon[1];
            U[INDEX(i, j, B2)] = bcon[2]*ucon[0] - bcon[0]*ucon[2];
            U[INDEX(i, j, B3)] = bcon[3]*ucon[0] - bcon[0]*ucon[3];

            U[INDEX(i, j, UU)] = tmp1*ucon[0]*ucov[0] + tmp2 - bcon[0]*bcov[0];
            U[INDEX(i, j, U1)] = tmp1*ucon[0]*ucov[1] - bcon[0]*bcov[1];
            U[INDEX(i, j, U2)] = tmp1*ucon[0]*ucov[2] - bcon[0]*bcov[2];
            U[INDEX(i, j, U3)] = tmp1*ucon[0]*ucov[3] - bcon[0]*bcov[3];
        }
    }

#undef ALPHA
#undef GCOV
#undef DELTA
}

void ComputedU_dt(const DOUBLE_SIZE* restrict prim,
                  const DOUBLE_SIZE* restrict dprim_dt,
                  const int iStart, const int iEnd,
                  const int jStart, const int jEnd,
                  DOUBLE_SIZE* restrict dU_dt)
{

#define VAR(varNum) (prim[INDEX(i, j, varNum)])
#define DVAR_DT(varNum) (dprim_dt[INDEX(i, j, varNum)])
#define SET_DU_DT(val, varNum) (dU_dt[varNum + DOF*(i+1)] = val)
#define DELTA(p, q) (p==q ? 1. : 0.)
#define GCOV(mu, nu) (mu==0 ? -DELTA(mu, nu) : DELTA(mu, nu))
#define ALPHA (1.)

    for (int j=jStart; j<jEnd; j++) {
        for (int i=iStart; i<iEnd; i++) {
            DOUBLE_SIZE gamma = sqrt(1 + GCOV(1, 1)*VAR(U1)*VAR(U1) + 
                                         GCOV(2, 2)*VAR(U2)*VAR(U2) +
                                         GCOV(3, 3)*VAR(U3)*VAR(U3) +
                                      2*(GCOV(1, 2)*VAR(U1)*VAR(U2) +
                                         GCOV(1, 3)*VAR(U1)*VAR(U3) +
                                         GCOV(2, 3)*VAR(U2)*VAR(U3)));

            DOUBLE_SIZE dgamma_dt = (2.*(GCOV(1, 1)*VAR(U1)*DVAR_DT(U1) + 
                                         GCOV(2, 2)*VAR(U2)*DVAR_DT(U2) +
                                         GCOV(3, 3)*VAR(U3)*DVAR_DT(U3)) +
                                     2.*(GCOV(1, 2)*DVAR_DT(U1)*VAR(U2) +
                                         GCOV(1, 2)*VAR(U1)*DVAR_DT(U2) +
                                         GCOV(1, 3)*DVAR_DT(U1)*VAR(U3) +
                                         GCOV(1, 3)*VAR(U1)*DVAR_DT(U3) +
                                         GCOV(2, 3)*DVAR_DT(U2)*VAR(U3) +
                                         GCOV(2, 3)*VAR(U2)*DVAR_DT(U3)))/gamma;

            DOUBLE_SIZE ucon[NDIM], ducon_dt[NDIM], ucov[NDIM], ducov_dt[NDIM];
            DOUBLE_SIZE bcon[NDIM], dbcon_dt[NDIM], bcov[NDIM], dbcov_dt[NDIM];

            ucon[0] = gamma/ALPHA;
            ucon[1] = VAR(U1) - gamma*GCOV(0, 1)*ALPHA;
            ucon[2] = VAR(U2) - gamma*GCOV(0, 2)*ALPHA;
            ucon[3] = VAR(U3) - gamma*GCOV(0, 3)*ALPHA;

            ducon_dt[0] = dgamma_dt/ALPHA;
            ducon_dt[1] = DVAR_DT(U1) - dgamma_dt*GCOV(0, 1)*ALPHA;
            ducon_dt[2] = DVAR_DT(U2) - dgamma_dt*GCOV(0, 2)*ALPHA;
            ducon_dt[3] = DVAR_DT(U3) - dgamma_dt*GCOV(0, 3)*ALPHA;

            ucov[0] = GCOV(0, 0)*ucon[0] +
                      GCOV(0, 1)*ucon[1] + 
                      GCOV(0, 2)*ucon[2] + 
                      GCOV(0, 3)*ucon[3];

            ucov[1] = GCOV(1, 0)*ucon[0] + 
                      GCOV(1, 1)*ucon[1] + 
                      GCOV(1, 2)*ucon[2] + 
                      GCOV(1, 3)*ucon[3];

            ucov[2] = GCOV(2, 0)*ucon[0] + 
                      GCOV(2, 1)*ucon[1] + 
                      GCOV(2, 2)*ucon[2] + 
                      GCOV(2, 3)*ucon[3];

            ucov[3] = GCOV(3, 0)*ucon[0] + 
                      GCOV(3, 1)*ucon[1] + 
                      GCOV(3, 2)*ucon[2] + 
                      GCOV(3, 3)*ucon[3];

            ducov_dt[0] = GCOV(0, 0)*ducon_dt[0] +
                          GCOV(0, 1)*ducon_dt[1] + 
                          GCOV(0, 2)*ducon_dt[2] + 
                          GCOV(0, 3)*ducon_dt[3];

            ducov_dt[1] = GCOV(1, 0)*ducon_dt[0] + 
                          GCOV(1, 1)*ducon_dt[1] + 
                          GCOV(1, 2)*ducon_dt[2] + 
                          GCOV(1, 3)*ducon_dt[3];

            ducov_dt[2] = GCOV(2, 0)*ducon_dt[0] + 
                          GCOV(2, 1)*ducon_dt[1] + 
                          GCOV(2, 2)*ducon_dt[2] + 
                          GCOV(2, 3)*ducon_dt[3];

            ducov_dt[3] = GCOV(3, 0)*ducon_dt[0] + 
                          GCOV(3, 1)*ducon_dt[1] + 
                          GCOV(3, 2)*ducon_dt[2] + 
                          GCOV(3, 3)*ducon_dt[3];


            bcon[0] = VAR(B1)*ucon[1] + VAR(B2)*ucon[2] + VAR(B3)*ucon[3];
            bcon[1] = (VAR(B1) + bcon[0]*ucon[1])/ucon[0];
            bcon[2] = (VAR(B2) + bcon[0]*ucon[2])/ucon[0];
            bcon[3] = (VAR(B3) + bcon[0]*ucon[3])/ucon[0];


            dbcon_dt[0] = DVAR_DT(B1)*ucon[1] + DVAR_DT(B2)*ucon[2] + 
                          DVAR_DT(B3)*ucon[3] + VAR(B1)*ducon_dt[1] + 
                          VAR(B2)*ducon_dt[2] + VAR(B3)*ducon_dt[3];

            dbcon_dt[1] = -(VAR(B1) + bcon[0]*ucon[1])*ducon_dt[0]/
                           (ucon[0]*ucon[0]) + 
                           (DVAR_DT(B1) + bcon[0]*ducon_dt[1] + 
                            dbcon_dt[0]*ucon[1])/ucon[0];
    
            dbcon_dt[2] = -(VAR(B2) + bcon[0]*ucon[2])*ducon_dt[0]/
                           (ucon[0]*ucon[0]) + 
                           (DVAR_DT(B2) + bcon[0]*ducon_dt[2] + 
                            dbcon_dt[0]*ucon[2])/ucon[0];

            dbcon_dt[3] = -(VAR(B3) + bcon[0]*ucon[3])*ducon_dt[0]/
                           (ucon[0]*ucon[0]) + 
                           (DVAR_DT(B3) + bcon[0]*ducon_dt[3] + 
                            dbcon_dt[0]*ucon[3])/ucon[0];

            bcov[0] = GCOV(0, 0)*bcon[0] +
                      GCOV(0, 1)*bcon[1] + 
                      GCOV(0, 2)*bcon[2] + 
                      GCOV(0, 3)*bcon[3];

            bcov[1] = GCOV(1, 0)*bcon[0] + 
                      GCOV(1, 1)*bcon[1] +
                      GCOV(1, 2)*bcon[2] + 
                      GCOV(1, 3)*bcon[3];

            bcov[2] = GCOV(2, 0)*bcon[0] + 
                      GCOV(2, 1)*bcon[1] + 
                      GCOV(2, 2)*bcon[2] + 
                      GCOV(2, 3)*bcon[3];

            bcov[3] = GCOV(3, 0)*bcon[0] + 
                      GCOV(3, 1)*bcon[1] + 
                      GCOV(3, 2)*bcon[2] + 
                      GCOV(3, 3)*bcon[3];

            dbcov_dt[0] = GCOV(0, 0)*dbcon_dt[0] +
                          GCOV(0, 1)*dbcon_dt[1] + 
                          GCOV(0, 2)*dbcon_dt[2] + 
                          GCOV(0, 3)*dbcon_dt[3];

            dbcov_dt[1] = GCOV(1, 0)*dbcon_dt[0] + 
                          GCOV(1, 1)*dbcon_dt[1] + 
                          GCOV(1, 2)*dbcon_dt[2] + 
                          GCOV(1, 3)*dbcon_dt[3];

            dbcov_dt[2] = GCOV(2, 0)*dbcon_dt[0] + 
                          GCOV(2, 1)*dbcon_dt[1] + 
                          GCOV(2, 2)*dbcon_dt[2] + 
                          GCOV(2, 3)*dbcon_dt[3];

            dbcov_dt[3] = GCOV(3, 0)*dbcon_dt[0] + 
                          GCOV(3, 1)*dbcon_dt[1] + 
                          GCOV(3, 2)*dbcon_dt[2] + 
                          GCOV(3, 3)*dbcon_dt[3];


            DOUBLE_SIZE bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] + 
                               bcon[2]*bcov[2] + bcon[3]*bcov[3];

            DOUBLE_SIZE dbsqr_dt = bcon[0]*dbcov_dt[0] + dbcon_dt[0]*bcov[0] +
                                   bcon[1]*dbcov_dt[1] + dbcon_dt[1]*bcov[1] +
                                   bcon[2]*dbcov_dt[2] + dbcon_dt[2]*bcov[2] +
                                   bcon[3]*dbcov_dt[3] + dbcon_dt[3]*bcov[3];

            int dir = 0;

            SET_DU_DT(DVAR_DT(RHO)*ucon[dir] + VAR(RHO)*ducon_dt[dir], RHO);

            SET_DU_DT(DVAR_DT(B1), B1);

            SET_DU_DT(DVAR_DT(B2), B2);

            SET_DU_DT(DVAR_DT(B3), B3);

            DOUBLE_SIZE P = (GAMMA - 1.)*VAR(UU);
            DOUBLE_SIZE dP_dt = (GAMMA - 1.)*DVAR_DT(UU);

            DOUBLE_SIZE tmp1 = dP_dt + DVAR_DT(RHO) + DVAR_DT(UU) + dbsqr_dt;
            DOUBLE_SIZE tmp2 = P + VAR(RHO) + VAR(UU) + bsqr;

            SET_DU_DT(tmp1*ucon[dir]*ucov[0] + 
                      tmp2*(ducon_dt[dir]*ucov[0] + ucon[dir]*ducov_dt[0]) + 
                      dP_dt + 0.5*dbsqr_dt - 
                      dbcon_dt[dir]*bcov[0] - bcon[dir]*dbcov_dt[0], UU);

            SET_DU_DT(tmp1*ucon[dir]*ucov[1] + 
                      tmp2*(ducon_dt[dir]*ucov[1] + ucon[dir]*ducov_dt[1]) - 
                      dbcon_dt[dir]*bcov[1] - bcon[dir]*dbcov_dt[1], U1);
        
            SET_DU_DT(tmp1*ucon[dir]*ucov[2] + 
                      tmp2*(ducon_dt[dir]*ucov[2] + ucon[dir]*ducov_dt[2]) - 
                      dbcon_dt[dir]*bcov[2] - bcon[dir]*dbcov_dt[2], U2);

            SET_DU_DT(tmp1*ucon[dir]*ucov[3] + 
                      tmp2*(ducon_dt[dir]*ucov[3] + ucon[dir]*ducov_dt[3]) - 
                      dbcon_dt[dir]*bcov[3] - bcon[dir]*dbcov_dt[3], U3);

        }
    }

#undef VAR
#undef DVAR_DT
#undef SET_DU_DT
#undef ALPHA
}

#define TOTAL_X1_SIZE 512
#define TOTAL_X2_SIZE 512
__kernel void ResFunctionCL(__global const DOUBLE_SIZE* restrict rho, 
                            __global const DOUBLE_SIZE* restrict u,
                            __global const DOUBLE_SIZE* restrict u1,
                            __global const DOUBLE_SIZE* restrict u2,
                            __global const DOUBLE_SIZE* restrict u3,
                            __global const DOUBLE_SIZE* restrict b1,
                            __global const DOUBLE_SIZE* restrict b2,
                            __global const DOUBLE_SIZE* restrict b3,
                            __global const DOUBLE_SIZE* restrict drho_dt,
                            __global const DOUBLE_SIZE* restrict du_dt,
                            __global const DOUBLE_SIZE* restrict du1_dt,
                            __global const DOUBLE_SIZE* restrict du2_dt,
                            __global const DOUBLE_SIZE* restrict du3_dt,
                            __global const DOUBLE_SIZE* restrict db1_dt,
                            __global const DOUBLE_SIZE* restrict db2_dt,
                            __global const DOUBLE_SIZE* restrict db3_dt,
                            __global DOUBLE_SIZE* restrict F_rho,
                            __global DOUBLE_SIZE* restrict F_u,
                            __global DOUBLE_SIZE* restrict F_u1,
                            __global DOUBLE_SIZE* restrict F_u2,
                            __global DOUBLE_SIZE* restrict F_u3,
                            __global DOUBLE_SIZE* restrict F_b1,
                            __global DOUBLE_SIZE* restrict F_b2,
                            __global DOUBLE_SIZE* restrict F_b3)
{
    int rowsPerWorkItem = TOTAL_X2_SIZE/TOTAL_NUM_WORK_ITEMS;
    
    int itemID = get_global_id(0);

    DOUBLE_SIZE prim[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE dprim_dt[(1 + 2*NG)*(TILE_SIZE+2)*DOF];

    DOUBLE_SIZE dU_dt[(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE prim_lx[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE prim_rx[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE prim_ly[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE prim_ry[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE flux_lx[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE flux_rx[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE flux_ly[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE flux_ry[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE U_lx[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE U_rx[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE U_ly[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE U_ry[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE flux_x[(1 + 2*NG)*(TILE_SIZE+2)*DOF];
    DOUBLE_SIZE flux_y[(1 + 2*NG)*(TILE_SIZE+2)*DOF];

    int X2End = rowsPerWorkItem;

    if (itemID==TOTAL_NUM_WORK_ITEMS-1) X2End = rowsPerWorkItem - NG;

    int i, iTile, globali, globalj;
    int j, loadj, workj, iStart, iEnd, jStart, jEnd;

#define COPY(i, j, prim, globali, globalj, var, varNum) { \
    prim[INDEX(i, j, varNum)] = var[globali+globalj*TOTAL_X1_SIZE/VECTOR_SIZE];\
}

    for (iTile=1; iTile<TOTAL_X1_SIZE/VECTOR_SIZE-1; iTile+=TILE_SIZE) {
        
        for (j=0; j<2*NG; j++) {
            for (i=-1; i<TILE_SIZE+1; i++) {
                /*i=-1 and i=TILE_SIZE contain the ghost zones*/
                
                int globali = i + iTile;
                int globalj = j + itemID*rowsPerWorkItem;

                COPY(i, j, prim, globali, globalj, rho, RHO)
                COPY(i, j, prim, globali, globalj, u, UU)
                COPY(i, j, prim, globali, globalj, u1, U1)
                COPY(i, j, prim, globali, globalj, u2, U2)
                COPY(i, j, prim, globali, globalj, u3, U3)
                COPY(i, j, prim, globali, globalj, b1, B1)
                COPY(i, j, prim, globali, globalj, b2, B2)
                COPY(i, j, prim, globali, globalj, b3, B3)

                COPY(i, j, dprim_dt, globali, globalj, drho_dt, RHO)
                COPY(i, j, dprim_dt, globali, globalj, du_dt, UU)
                COPY(i, j, dprim_dt, globali, globalj, du1_dt, U1)
                COPY(i, j, dprim_dt, globali, globalj, du2_dt, U2)
                COPY(i, j, dprim_dt, globali, globalj, du3_dt, U3)
                COPY(i, j, dprim_dt, globali, globalj, db1_dt, B1)
                COPY(i, j, dprim_dt, globali, globalj, db2_dt, B2)
                COPY(i, j, dprim_dt, globali, globalj, db3_dt, B3)
            }
        }

        for (j=NG; j<X2End; j++) {

            loadj = j+NG;
            workj = j;

            for (i=-1; i<TILE_SIZE+1; i++) {

                globali = i + iTile;
                globalj = loadj + itemID*rowsPerWorkItem;

                COPY(i, loadj, prim, globali, globalj, rho, RHO)
                COPY(i, loadj, prim, globali, globalj, u, UU)
                COPY(i, loadj, prim, globali, globalj, u1, U1)
                COPY(i, loadj, prim, globali, globalj, u2, U2)
                COPY(i, loadj, prim, globali, globalj, u3, U3)
                COPY(i, loadj, prim, globali, globalj, b1, B1)
                COPY(i, loadj, prim, globali, globalj, b2, B2)
                COPY(i, loadj, prim, globali, globalj, b3, B3)

                COPY(i, loadj, dprim_dt, globali, globalj, drho_dt, RHO)
                COPY(i, loadj, dprim_dt, globali, globalj, du_dt, UU)
                COPY(i, loadj, dprim_dt, globali, globalj, du1_dt, U1)
                COPY(i, loadj, dprim_dt, globali, globalj, du2_dt, U2)
                COPY(i, loadj, dprim_dt, globali, globalj, du3_dt, U3)
                COPY(i, loadj, dprim_dt, globali, globalj, db1_dt, B1)
                COPY(i, loadj, dprim_dt, globali, globalj, db2_dt, B2)
                COPY(i, loadj, dprim_dt, globali, globalj, db3_dt, B3)

            }

            ComputedU_dt(prim, dprim_dt, -1, TILE_SIZE+1, workj, workj+1, dU_dt);

//            Reconstruct(prim, 
//                        -1, TILE_SIZE+1,
//                        workj-NG, workj+NG,
//                        prim_lx, prim_rx,
//                        prim_ly, prim_ry);

            ComputeFluxAndU(prim_lx, 1, -1, TILE_SIZE+1,
                            workj-NG, workj+NG, flux_lx, U_lx);

            ComputeFluxAndU(prim_rx, 1, -1, TILE_SIZE+1,
                            workj-NG, workj+NG, flux_rx, U_rx);
            
            ComputeFluxAndU(prim_ly, 2, -1, TILE_SIZE+1, 
                            workj-NG, workj+NG, flux_ly, U_ly);

            ComputeFluxAndU(prim_ry, 2, -1, TILE_SIZE+1, 
                            workj-NG, workj+NG, flux_ry, U_ry);

//            RiemannSolver(flux_lx, flux_rx,
//                          flux_ly, flux_ry,
//                          U_lx, U_rx,
//                          U_ly, U_ry,
//                          -1, TILE_SIZE+1,
//                          workj-NG, workj+NG,
//                          flux_x, flux_y);
//
//#define COMPUTE_F(i, j, iOffset, globali, globalj, var, varNum) { \
//    var[globali + globalj*TOTAL_X1_SIZE/VECTOR_SIZE] = \
//        dU_dt[varNum + DOF*(i+1)] + \
//        (LOAD_VAR(i, j, iOffset+1, 0, flux_x, varNum) - \
//         LOAD_VAR(i, j, iOffset, 0, flux_x, varNum))/DX1 + \
//        (LOAD_VAR(i, j, iOffset, 1, flux_y, varNum) - \
//         LOAD_VAR(i, j, iOffset, 0, flux_y, varNum))/DX2; \
//}
#define COMPUTE_F(i, j, iOffset, globali, globalj, var, varNum) { \
    var[globali + globalj*TOTAL_X1_SIZE/VECTOR_SIZE] = \
        dU_dt[varNum + DOF*(i+1)]; \
}

            #pragma vector nontemporal
            for (i=-1; i<TILE_SIZE; i++) {
                globali = i + iTile;
                globalj = workj + itemID*rowsPerWorkItem;
                
                COMPUTE_F(i, workj, 0, globali, globalj, F_rho, RHO);
                COMPUTE_F(i, workj, 0, globali, globalj, F_u, UU);
                COMPUTE_F(i, workj, 0, globali, globalj, F_u1, U1);
                COMPUTE_F(i, workj, 0, globali, globalj, F_u2, U2);
                COMPUTE_F(i, workj, 0, globali, globalj, F_u3, U3);
                COMPUTE_F(i, workj, 0, globali, globalj, F_b1, B1);
                COMPUTE_F(i, workj, 0, globali, globalj, F_b2, B2);
                COMPUTE_F(i, workj, 0, globali, globalj, F_b3, B3);
            }

        }

    }

#undef COPY
#undef COMPUTE_F
}


