#define ADIABATIC_INDEX 4./3.

#define DELTA(p, q) (p==q ? 1. : 0.)
#define GCOV(mu, nu) (mu==0 ? -DELTA(mu, nu) : DELTA(mu, nu))
#define ALPHA (1.)

#define GAMMA(var) (sqrt(1 + GCOV(1, 1)*var[U1]*var[U1] + \
                             GCOV(2, 2)*var[U2]*var[U2] + \
                             GCOV(3, 3)*var[U3]*var[U3] + \
                          2*(GCOV(1, 2)*var[U1]*var[U2] + \
                             GCOV(1, 3)*var[U1]*var[U3] + \
                             GCOV(2, 3)*var[U2]*var[U3])) )

#define DGAMMA_DT(var,dvar_dt) ( (2.*(GCOV(1, 1)*var[U1]*dvar_dt[U1] + \
                                      GCOV(2, 2)*var[U2]*dvar_dt[U2] + \
                                      GCOV(3, 3)*var[U3]*dvar_dt[U3]) + \
                                  2.*(GCOV(1, 2)*dvar_dt[U1]*var[U2] + \
                                      GCOV(1, 2)*var[U1]*dvar_dt[U2] + \
                                      GCOV(1, 3)*dvar_dt[U1]*var[U3] + \
                                      GCOV(1, 3)*var[U1]*dvar_dt[U3] + \
                                      GCOV(2, 3)*dvar_dt[U2]*var[U3] + \
                                      GCOV(2, 3)*var[U2]*dvar_dt[U3]))/ \
                                      GAMMA(var))

#define UCON(var,mu) (mu==0 ? GAMMA(var)/ALPHA : \
                      var[UU+mu] - GAMMA(var)*GCOV(0,mu)*ALPHA)

#define DUCON_DT(var,dvar_dt,mu) (mu==0 ? DGAMMA_DT(var,dvar_dt)/ALPHA : \
                                  dvar_dt[UU+mu] - \
                                  DGAMMA_DT(var,dvar_dt)*GCOV(0, mu)*ALPHA)

#define UCOV(var,mu) (GCOV(mu, 0)*UCON(var,0) + GCOV(mu, 1)*UCON(var,1) + \
                      GCOV(mu, 2)*UCON(var,2) + GCOV(mu, 3)*UCON(var,3) )

#define DUCOV_DT(var,dvar_dt,mu) (GCOV(mu, 0)*DUCON_DT(var,dvar_dt,0) + \
                                  GCOV(mu, 1)*DUCON_DT(var,dvar_dt,1) + \
                                  GCOV(mu, 2)*DUCON_DT(var,dvar_dt,2) + \
                                  GCOV(mu, 3)*DUCON_DT(var,dvar_dt,3) )

#define BCON0(var) (var[B1]*UCON(var,1)+var[B2]*UCON(var,2)+var[B3]*UCON(var,3))

#define BCON(var,mu) (mu==0 ? BCON0(var) : \
                      (var[U3+mu] + BCON0(var)*UCON(var,mu) )/UCON(var,0) )

#define DBCON_DT0(var,dvar_dt) (dvar_dt[B1]*UCON(var,1) + \
                                dvar_dt[B2]*UCON(var,2) + \
                                dvar_dt[B3]*UCON(var,3) + \
                                var[B1]*DUCON_DT(var,dvar_dt,1) + \
                                var[B3]*DUCON_DT(var,dvar_dt,2) + \
                                var[B2]*DUCON_DT(var,dvar_dt,3) )

#define DBCON_DT(var,dvar_dt,mu) (mu==0 ? DBCON_DT0(var,dvar_dt) : \
                                  -(var[U3+mu] + \
                                    BCON0(var)*UCON(var,mu))*\
                                  DUCON_DT(var,dvar_dt,0)/ \
                                  (UCON(var,0)*UCON(var,0)) + \
                                  (dvar_dt[U3+mu] + BCON0(var)*\
                                   DUCON_DT(var,dvar_dt,mu) + \
                                   DBCON_DT0(var,dvar_dt)*UCON(var,mu))/\
                                  UCON(var,0))

#define BCOV(var,mu) (GCOV(mu, 0)*BCON(var,0) + GCOV(mu, 1)*BCON(var,1) + \
                      GCOV(mu, 2)*BCON(var,2) + GCOV(mu, 3)*BCON(var,3) )

#define DBCOV_DT(var,dvar_dt,mu) (GCOV(mu, 0)*DBCON_DT(var,dvar_dt,0) + \
                                  GCOV(mu, 1)*DBCON_DT(var,dvar_dt,1) + \
                                  GCOV(mu, 2)*DBCON_DT(var,dvar_dt,2) + \
                                  GCOV(mu, 3)*DBCON_DT(var,dvar_dt,3) )

#define BSQR(var) (BCON(var,0)*BCOV(var,0) + BCON(var,1)*BCOV(var,1) + \
                   BCON(var,2)*BCOV(var,2) + BCON(var,3)*BCOV(var,3) )

#define DBSQR_DT(var,dvar_dt) (BCON(var,0)*DBCOV_DT(var,dvar_dt,0) + \
                               DBCON_DT(var,dvar_dt,0)*BCOV(var,0) + \
                               BCON(var,1)*DBCOV_DT(var,dvar_dt,1) + \
                               DBCON_DT(var,dvar_dt,1)*BCOV(var,1) + \
                               BCON(var,2)*DBCOV_DT(var,dvar_dt,2) + \
                               DBCON_DT(var,dvar_dt,2)*BCOV(var,2) + \
                               BCON(var,3)*DBCOV_DT(var,dvar_dt,3) + \
                               DBCON_DT(var,dvar_dt,3)*BCOV(var,3) )

#define P(var) ((ADIABATIC_INDEX - 1.)*var[UU])
#define DP_DT(var,dvar_dt) ((ADIABATIC_INDEX - 1.)*dvar_dt[UU])

#define TMP1(var) (P(var) + var[RHO] + var[UU] + BSQR(var))
#define DTMP1_DT(var,dvar_dt) (DP_DT(var,dvar_dt) + \
                               dvar_dt[RHO] + dvar_dt[UU] + \
                               DBSQR_DT(var,dvar_dt))

#define TMP2(var) (P(var) + 0.5*BSQR(var))

//void ComputeFluxAndU(const REAL* restrict var,
//                     REAL* restrict flux,
//                     REAL* restrict U,
//                     const int dir) 
//{
//    flux[RHO]  = var[RHO]*UCON(var, dir);
//    flux[UU] = TMP1(var)*UCON(var,dir)*UCOV(var,0) +
//               TMP2(var)*DELTA(dir, 0)
//               - BCON(var,dir)*BCOV(var,0);
//    flux[U1] = TMP1(var)*UCON(var,dir)*UCOV(var,1) +
//               TMP2(var)*DELTA(dir, 1)
//               - BCON(var,dir)*BCOV(var,1);
//    flux[U2] = TMP1(var)*UCON(var,dir)*UCOV(var,2) +
//               TMP2(var)*DELTA(dir, 2)
//               - BCON(var,dir)*BCOV(var,2);
//    flux[U3] = TMP1(var)*UCON(var,dir)*UCOV(var,3) +
//               TMP2(var)*DELTA(dir, 3)
//               - BCON(var,dir)*BCOV(var,3);
//    flux[B1] = BCON(var,1)*UCON(var,dir) - BCON(var,dir)*UCON(var,1);
//    flux[B2] = BCON(var,2)*UCON(var,dir) - BCON(var,dir)*UCON(var,2);
//    flux[B3] = BCON(var,3)*UCON(var,dir) - BCON(var,dir)*UCON(var,3);
//    U[RHO]  = var[RHO]*UCON(var,0);
//    U[B1] = BCON(var,1)*UCON(var,0) - BCON(var,0)*UCON(var,1);
//    U[B2] = BCON(var,2)*UCON(var,0) - BCON(var,0)*UCON(var,2);
//    U[B3] = BCON(var,3)*UCON(var,0) - BCON(var,0)*UCON(var,3);
//    U[UU] = TMP1(var)*UCON(var,0)*UCOV(var,0) +
//            TMP2(var) - BCON(var,0)*BCOV(var,0);
//    U[U1] = TMP1(var)*UCON(var,0)*UCOV(var,1) - BCON(var,0)*BCOV(var,1);
//    U[U2] = TMP1(var)*UCON(var,0)*UCOV(var,2) - BCON(var,0)*BCOV(var,2);
//    U[U3] = TMP1(var)*UCON(var,0)*UCOV(var,3) - BCON(var,0)*BCOV(var,3);
//}
//
//void ComputedU_dt(const REAL* restrict var,
//                  const REAL* restrict dvar_dt,
//                  REAL* restrict dU_dt)
//{
//    dU_dt[RHO] = dvar_dt[RHO]*UCON(var,0) + var[RHO]*DUCON_DT(var,dvar_dt,0);
//    dU_dt[B1] = dvar_dt[B1];
//    dU_dt[B2] = dvar_dt[B2];
//    dU_dt[B3] = dvar_dt[B3];
//    dU_dt[UU] = DTMP1_DT(var,dvar_dt)*UCON(var,0)*UCOV(var,0) +
//                TMP1(var)*(DUCON_DT(var,dvar_dt,0)*UCOV(var,0) +
//                UCON(var,0)*DUCOV_DT(var,dvar_dt,0)) +
//                DP_DT(var,dvar_dt) + 0.5*DBSQR_DT(var,dvar_dt) -
//                DBCON_DT(var,dvar_dt,0)*BCOV(var,0) -
//                BCON(var,0)*DBCOV_DT(var,dvar_dt,0);
//    dU_dt[U1] = DTMP1_DT(var,dvar_dt)*UCON(var,0)*UCOV(var,1) +
//                TMP1(var)*(DUCON_DT(var,dvar_dt,0)*UCOV(var,1) +
//                UCON(var,0)*DUCOV_DT(var,dvar_dt,1)) -
//                DBCON_DT(var,dvar_dt,0)*BCOV(var,1) -
//                BCON(var,0)*DBCOV_DT(var,dvar_dt,1);
//    dU_dt[U2] = DTMP1_DT(var,dvar_dt)*UCON(var,0)*UCOV(var,2) +
//                TMP1(var)*(DUCON_DT(var,dvar_dt,0)*UCOV(var,2) +
//                UCON(var,0)*DUCOV_DT(var,dvar_dt,2)) -
//                DBCON_DT(var,dvar_dt,0)*BCOV(var,2) -
//                BCON(var,0)*DBCOV_DT(var,dvar_dt,2);
//    dU_dt[U3] = DTMP1_DT(var,dvar_dt)*UCON(var,0)*UCOV(var,3) +
//                TMP1(var)*(DUCON_DT(var,dvar_dt,0)*UCOV(var,3) +
//                UCON(var,0)*DUCOV_DT(var,dvar_dt,3)) -
//                DBCON_DT(var,dvar_dt,0)*BCOV(var,3) -
//                BCON(var,0)*DBCOV_DT(var,dvar_dt,3);
//}


void ComputeFluxAndU(const REAL* restrict var,
                     REAL* restrict flux,
                     REAL* restrict U,
                     const int dir) 
{
    REAL gamma, ucon[NDIM], ucov[NDIM], bcov[NDIM], bcon[NDIM];
    REAL tmp1, tmp2, bsqr;

    gamma = sqrt(1 + GCOV(1, 1)*var[U1]*var[U1] +
                     GCOV(2, 2)*var[U2]*var[U2] +
                     GCOV(3, 3)*var[U3]*var[U3] +
                  2*(GCOV(1, 2)*var[U1]*var[U2] +
                     GCOV(1, 3)*var[U1]*var[U3] +
                     GCOV(2, 3)*var[U2]*var[U3]));

    ucon[0] = gamma/ALPHA;
    ucon[1] = var[U1] - gamma*GCOV(0,1)*ALPHA;
    ucon[2] = var[U2] - gamma*GCOV(0,2)*ALPHA;
    ucon[3] = var[U3] - gamma*GCOV(0,3)*ALPHA;

    ucov[0] = GCOV(0, 0)*ucon[0] + GCOV(0, 1)*ucon[1] +
              GCOV(0, 2)*ucon[2] + GCOV(0, 3)*ucon[3];
    ucov[1] = GCOV(1, 0)*ucon[0] + GCOV(1, 1)*ucon[1] +
              GCOV(1, 2)*ucon[2] + GCOV(1, 3)*ucon[3];
    ucov[2] = GCOV(2, 0)*ucon[0] + GCOV(2, 1)*ucon[1] +
              GCOV(2, 2)*ucon[2] + GCOV(2, 3)*ucon[3];
    ucov[3] = GCOV(3, 0)*ucon[0] + GCOV(3, 1)*ucon[1] +
              GCOV(3, 2)*ucon[2] + GCOV(3, 3)*ucon[3];

    bcon[0] = var[B1]*ucon[1] + var[B2]*ucon[2] + var[B3]*ucon[3];
    bcon[1] = (var[B1] + bcon[0]*ucon[1])/ucon[0];
    bcon[2] = (var[B2] + bcon[0]*ucon[2])/ucon[0];
    bcon[3] = (var[B3] + bcon[0]*ucon[3])/ucon[0];

    bcov[0] = GCOV(0, 0)*bcon[0] + GCOV(0, 1)*bcon[1] +
              GCOV(0, 2)*bcon[2] + GCOV(0, 3)*bcon[3];
    bcov[1] = GCOV(1, 0)*bcon[0] + GCOV(1, 1)*bcon[1] +
              GCOV(1, 2)*bcon[2] + GCOV(1, 3)*bcon[3];
    bcov[2] = GCOV(2, 0)*bcon[0] + GCOV(2, 1)*bcon[1] +
              GCOV(2, 2)*bcon[2] + GCOV(2, 3)*bcon[3];
    bcov[3] = GCOV(3, 0)*bcon[0] + GCOV(3, 1)*bcon[1] +
              GCOV(3, 2)*bcon[2] + GCOV(3, 3)*bcon[3];

    bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] + 
           bcon[2]*bcov[2] + bcon[3]*bcov[3];

    tmp1 = P(var) + var[RHO] + var[UU] + bsqr;
    tmp2 = P(var) + 0.5*bsqr;

    flux[RHO]  = var[RHO]*ucon[dir];

    flux[UU] = tmp1*ucon[dir]*ucov[0] + tmp2*DELTA(dir, 0)
               - bcon[dir]*bcov[0];

    flux[U1] = tmp1*ucon[dir]*ucov[1] + tmp2*DELTA(dir, 1)
               - bcon[dir]*bcov[1];

    flux[U2] = tmp1*ucon[dir]*ucov[2] + tmp2*DELTA(dir, 2)
               - bcon[dir]*bcov[2];

    flux[U3] = tmp1*ucon[dir]*ucov[3] + tmp2*DELTA(dir, 3)
               - bcon[dir]*bcov[3];

    flux[B1] = bcon[1]*ucon[dir] - bcon[dir]*ucon[1];
    flux[B2] = bcon[2]*ucon[dir] - bcon[dir]*ucon[2];
    flux[B3] = bcon[3]*ucon[dir] - bcon[dir]*ucon[3];

    U[RHO]  = var[RHO]*ucon[0];

    U[B1] = bcon[1]*ucon[0] - bcon[0]*ucon[1];
    U[B2] = bcon[2]*ucon[0] - bcon[0]*ucon[2];
    U[B3] = bcon[3]*ucon[0] - bcon[0]*ucon[3];

    U[UU] = tmp1*ucon[0]*ucov[0] + tmp2 - bcon[0]*bcov[0];
    U[U1] = tmp1*ucon[0]*ucov[1] - bcon[0]*bcov[1];
    U[U2] = tmp1*ucon[0]*ucov[2] - bcon[0]*bcov[2];
    U[U3] = tmp1*ucon[0]*ucov[3] - bcon[0]*bcov[3];
}


void ComputedU_dt(const REAL* restrict var,
                  const REAL* restrict dvar_dt,
                  REAL* restrict dU_dt)
{
    REAL gamma, dgamma_dt, tmp1, tmp2, bsqr, dbsqr_dt;
    REAL ucon[NDIM], ducon_dt[NDIM], ucov[NDIM], ducov_dt[NDIM];
    REAL bcon[NDIM], dbcon_dt[NDIM], bcov[NDIM], dbcov_dt[NDIM];

    gamma = sqrt(1 + GCOV(1, 1)*var[U1]*var[U1] + 
                     GCOV(2, 2)*var[U2]*var[U2] +
                     GCOV(3, 3)*var[U3]*var[U3] +
                  2*(GCOV(1, 2)*var[U1]*var[U2] +
                     GCOV(1, 3)*var[U1]*var[U3] +
                     GCOV(2, 3)*var[U2]*var[U3]));

    dgamma_dt = (2.*(GCOV(1, 1)*var[U1]*dvar_dt[U1] + 
                     GCOV(2, 2)*var[U2]*dvar_dt[U2] +
                     GCOV(3, 3)*var[U3]*dvar_dt[U3]) +
                 2.*(GCOV(1, 2)*dvar_dt[U1]*var[U2] +
                     GCOV(1, 2)*var[U1]*dvar_dt[U2] +
                     GCOV(1, 3)*dvar_dt[U1]*var[U3] +
                     GCOV(1, 3)*var[U1]*dvar_dt[U3] +
                     GCOV(2, 3)*dvar_dt[U2]*var[U3] +
                     GCOV(2, 3)*var[U2]*dvar_dt[U3]))/gamma;

    ucon[0] = gamma/ALPHA;
    ucon[1] = var[U1] - gamma*GCOV(0, 1)*ALPHA;
    ucon[2] = var[U2] - gamma*GCOV(0, 2)*ALPHA;
    ucon[3] = var[U3] - gamma*GCOV(0, 3)*ALPHA;

    ducon_dt[0] = dgamma_dt/ALPHA;
    ducon_dt[1] = dvar_dt[U1] - dgamma_dt*GCOV(0, 1)*ALPHA;
    ducon_dt[2] = dvar_dt[U2] - dgamma_dt*GCOV(0, 2)*ALPHA;
    ducon_dt[3] = dvar_dt[U3] - dgamma_dt*GCOV(0, 3)*ALPHA;

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


    bcon[0] = var[B1]*ucon[1] + var[B2]*ucon[2] + var[B3]*ucon[3];
    bcon[1] = (var[B1] + bcon[0]*ucon[1])/ucon[0];
    bcon[2] = (var[B2] + bcon[0]*ucon[2])/ucon[0];
    bcon[3] = (var[B3] + bcon[0]*ucon[3])/ucon[0];


    dbcon_dt[0] = dvar_dt[B1]*ucon[1] + dvar_dt[B2]*ucon[2] + 
                  dvar_dt[B3]*ucon[3] + var[B1]*ducon_dt[1] + 
                  var[B2]*ducon_dt[2] + var[B3]*ducon_dt[3];

    dbcon_dt[1] = -(var[B1] + bcon[0]*ucon[1])*ducon_dt[0]/
                   (ucon[0]*ucon[0]) + 
                   (dvar_dt[B1] + bcon[0]*ducon_dt[1] + 
                    dbcon_dt[0]*ucon[1])/ucon[0];
    
    dbcon_dt[2] = -(var[B2] + bcon[0]*ucon[2])*ducon_dt[0]/
                   (ucon[0]*ucon[0]) + 
                   (dvar_dt[B2] + bcon[0]*ducon_dt[2] + 
                    dbcon_dt[0]*ucon[2])/ucon[0];

    dbcon_dt[3] = -(var[B3] + bcon[0]*ucon[3])*ducon_dt[0]/
                   (ucon[0]*ucon[0]) + 
                   (dvar_dt[B3] + bcon[0]*ducon_dt[3] + 
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


    bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] + 
           bcon[2]*bcov[2] + bcon[3]*bcov[3];

    dbsqr_dt = bcon[0]*dbcov_dt[0] + dbcon_dt[0]*bcov[0] +
               bcon[1]*dbcov_dt[1] + dbcon_dt[1]*bcov[1] +
               bcon[2]*dbcov_dt[2] + dbcon_dt[2]*bcov[2] +
               bcon[3]*dbcov_dt[3] + dbcon_dt[3]*bcov[3];

    tmp1 = DP_DT(var, dvar_dt) + dvar_dt[RHO] + dvar_dt[UU] + dbsqr_dt;
    tmp2 = P(var) + var[RHO] + var[UU] + bsqr;

    dU_dt[RHO] = dvar_dt[RHO]*ucon[0] + var[RHO]*ducon_dt[0];

    dU_dt[B1] = dvar_dt[B1];
    dU_dt[B2] = dvar_dt[B2];
    dU_dt[B3] = dvar_dt[B3];

    dU_dt[UU] = tmp1*ucon[0]*ucov[0] + 
                tmp2*(ducon_dt[0]*ucov[0] + ucon[0]*ducov_dt[0]) + 
                DP_DT(var,dvar_dt) + 0.5*dbsqr_dt - 
                dbcon_dt[0]*bcov[0] - bcon[0]*dbcov_dt[0];

    dU_dt[U1] = tmp1*ucon[0]*ucov[1] + 
                tmp2*(ducon_dt[0]*ucov[1] + ucon[0]*ducov_dt[1]) - 
                dbcon_dt[0]*bcov[1] - bcon[0]*dbcov_dt[1];
    
    dU_dt[U2] = tmp1*ucon[0]*ucov[2] + 
                tmp2*(ducon_dt[0]*ucov[2] + ucon[0]*ducov_dt[2]) - 
                dbcon_dt[0]*bcov[2] - bcon[0]*dbcov_dt[2];

    dU_dt[U3] = tmp1*ucon[0]*ucov[3] + 
                tmp2*(ducon_dt[0]*ucov[3] + ucon[0]*ducov_dt[3]) - 
                dbcon_dt[0]*bcov[3] - bcon[0]*dbcov_dt[3];
}
