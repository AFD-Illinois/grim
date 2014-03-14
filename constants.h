#define COMPUTE_DIM 2
#define NDIM 4
#define N1 64
#define N2 64
#define NG 2

#define REAL double
#define TILE_SIZE_X1 8
#define TILE_SIZE_X2 8

#define A_SPIN 0.9375
#define M 1.
#define R_MIN 6.
#define R_MAX 12.
#define H_SLOPE 0.3
#define R0 0.
#define R_IN .98*(1. + sqrt(1. - A_SPIN*A_SPIN))
#define R_OUT 40.
#define X1_START log(R_IN - R0)
#define X2_START 0.1
#define DX1 (log((R_OUT - R0)/(R_IN - R0))/(REAL)N1)
#define DX2 ((1.-2.*X2_START)/(REAL)N2)
//#define DX1 (1./(REAL)N1)
//#define DX2 (1./(REAL)N2)
#define DT 1e-5
#define KAPPA 1e-3
#define BETA 1e2
#define ADIABATIC_INDEX (4./3.)
#define RHO_MIN (1e-5)
#define U_MIN (1e-7)
#define RHO_MIN_LIMIT (1e-15)
#define U_MIN_LIMIT (1e-15)
#define GAMMA_MAX (50.)

#define EPS (1e-5)

#define RHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7
#define DOF 8

// iTile, jTile have ranges [-NG, TILE_SIZE+NG)
#define INDEX_LOCAL(iTile,jTile,var) (iTile+NG + \
                                      (TILE_SIZE_X1+2*NG)*(jTile+NG + \
                                      (TILE_SIZE_X2+2*NG)*(var)))
//#define INDEX_LOCAL(iTile,jTile,var) ((var) + \
//                                      DOF*((iTile)+NG + \
//                                      (TILE_SIZE_X1+2*NG)*((jTile)+NG)))

// i, j have ranges [0, N1), [0, N2)
#define INDEX_GLOBAL(i,j,var) (var + DOF*((i)+(N1)*(j)))

#define i_TO_X1_CENTER(i) (X1_START + (i + 0.5)*DX1)
#define j_TO_X2_CENTER(j) (X2_START + (j + 0.5)*DX2)
#define i_TO_X1_FACE(i) (X1_START + (i)*DX1)
#define j_TO_X2_FACE(j) (X2_START + (j)*DX2)

void BLCoords(REAL* r, REAL* theta,
              const REAL X1, const REAL X2)
{
    *r = exp(X1) + R0;
    *theta = M_PI*(X2) + ((1 - H_SLOPE)/2.)*sin(2.*M_PI*(X2));
}


void gCovCalc(REAL gcov[NDIM][NDIM],
              const REAL X1, const REAL X2)
{
    REAL r, theta;
    BLCoords(&r, &theta, X1, X2);

    REAL rhoSqr = r*r + A_SPIN*A_SPIN*cos(theta)*cos(theta);
    REAL rFactor = r - R0;
    REAL hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*(X2));

    gcov[0][0] = -1. + 2.*r/rhoSqr;
    gcov[0][1] = (2.*r/rhoSqr)*rFactor;
    gcov[0][2] = 0.;
    gcov[0][3] = -2.*A_SPIN*r*sin(theta)*sin(theta)/rhoSqr;

    gcov[1][0] = gcov[0][1];
    gcov[1][1] = (1. + 2.*r/rhoSqr)*rFactor*rFactor;
    gcov[1][2] = 0.;
    gcov[1][3] = -A_SPIN*sin(theta)*sin(theta)*(1.+
                2.*r/rhoSqr)*rFactor;

    gcov[2][0] = gcov[0][2];
    gcov[2][1] = gcov[1][2];
    gcov[2][2] = rhoSqr*hFactor*hFactor;
    gcov[2][3] = 0.;


    gcov[3][0] = gcov[0][3];
    gcov[3][1] = gcov[1][3];
    gcov[3][2] = gcov[2][3];
    gcov[3][3] = sin(theta)*sin(theta)*(rhoSqr +
                            A_SPIN*A_SPIN*sin(theta)*sin(theta)*\
                            (1. + 2.*r/rhoSqr) );

//    gcov[0][0] = -1.;
//    gcov[0][1] = 0.;
//    gcov[0][2] = 0.;
//    gcov[0][3] = 0.;
//
//    gcov[1][0] = 0.;
//    gcov[1][1] = 1.;
//    gcov[1][2] = 0.;
//    gcov[1][3] = 0.;
//
//    gcov[2][0] = 0.;
//    gcov[2][1] = 0.;
//    gcov[2][2] = 1.;
//    gcov[2][3] = 0.;
//
//    gcov[3][0] = 0.;
//    gcov[3][1] = 0.;
//    gcov[3][2] = 0.;
//    gcov[3][3] = 1.;
}

void gConCalc(REAL gcon[NDIM][NDIM],
              const REAL gcov[NDIM][NDIM],
              const REAL gdet)
{

    gcon[0][0] = 
        (gcov[1][1]*gcov[2][2]*gcov[3][3] + 
         gcov[1][2]*gcov[2][3]*gcov[3][1] + 
         gcov[1][3]*gcov[2][1]*gcov[3][2] - 
         gcov[1][1]*gcov[2][3]*gcov[3][2] - 
         gcov[1][2]*gcov[2][1]*gcov[3][3] - 
         gcov[1][3]*gcov[2][2]*gcov[3][1])/gdet;

    gcon[0][1] = 
        (gcov[0][1]*gcov[2][3]*gcov[3][2] + 
         gcov[0][2]*gcov[2][1]*gcov[3][3] + 
         gcov[0][3]*gcov[2][2]*gcov[3][1] - 
         gcov[0][1]*gcov[2][2]*gcov[3][3] - 
         gcov[0][2]*gcov[2][3]*gcov[3][1] - 
         gcov[0][3]*gcov[2][1]*gcov[3][2])/gdet;

    gcon[0][2] = 
        (gcov[0][1]*gcov[1][2]*gcov[3][3] + 
         gcov[0][2]*gcov[1][3]*gcov[3][1] + 
         gcov[0][3]*gcov[1][1]*gcov[3][2] - 
         gcov[0][1]*gcov[1][3]*gcov[3][2] - 
         gcov[0][2]*gcov[1][1]*gcov[3][3] - 
         gcov[0][3]*gcov[1][2]*gcov[3][1])/gdet;

    gcon[0][3] = 
        (gcov[0][1]*gcov[1][3]*gcov[2][2] + 
         gcov[0][2]*gcov[1][1]*gcov[2][3] + 
         gcov[0][3]*gcov[1][2]*gcov[2][1] - 
         gcov[0][1]*gcov[1][2]*gcov[2][3] - 
         gcov[0][2]*gcov[1][3]*gcov[2][1] - 
         gcov[0][3]*gcov[1][1]*gcov[2][2])/gdet;

    gcon[1][0] = gcon[0][1];
    
    gcon[1][1] = 
        (gcov[0][0]*gcov[2][2]*gcov[3][3] + 
         gcov[0][2]*gcov[2][3]*gcov[3][0] + 
         gcov[0][3]*gcov[2][0]*gcov[3][2] - 
         gcov[0][0]*gcov[2][3]*gcov[3][2] - 
         gcov[0][2]*gcov[2][0]*gcov[3][3] - 
         gcov[0][3]*gcov[2][2]*gcov[3][0])/gdet;

    gcon[1][2] = 
        (gcov[0][0]*gcov[1][3]*gcov[3][2] + 
         gcov[0][2]*gcov[1][0]*gcov[3][3] + 
         gcov[0][3]*gcov[1][2]*gcov[3][0] - 
         gcov[0][0]*gcov[1][2]*gcov[3][3] - 
         gcov[0][2]*gcov[1][3]*gcov[3][0] - 
         gcov[0][3]*gcov[1][0]*gcov[3][2])/gdet;

    gcon[1][3] = 
        (gcov[0][0]*gcov[1][2]*gcov[2][3] + 
         gcov[0][2]*gcov[1][3]*gcov[2][0] + 
         gcov[0][3]*gcov[1][0]*gcov[2][2] - 
         gcov[0][0]*gcov[1][3]*gcov[2][2] - 
         gcov[0][2]*gcov[1][0]*gcov[2][3] - 
         gcov[0][3]*gcov[1][2]*gcov[2][0])/gdet;

    gcon[2][0] = gcon[0][2];
    gcon[2][1] = gcon[1][2];

    gcon[2][2] =
        (gcov[0][0]*gcov[1][1]*gcov[3][3] + 
         gcov[0][1]*gcov[1][3]*gcov[3][0] + 
         gcov[0][3]*gcov[1][0]*gcov[3][1] - 
         gcov[0][0]*gcov[1][3]*gcov[3][1] - 
         gcov[0][1]*gcov[1][0]*gcov[3][3] - 
         gcov[0][3]*gcov[1][1]*gcov[3][0])/gdet;

    gcon[2][3] =
        (gcov[0][0]*gcov[1][3]*gcov[2][1] + 
         gcov[0][1]*gcov[1][0]*gcov[2][3] + 
         gcov[0][3]*gcov[1][1]*gcov[2][0] - 
         gcov[0][0]*gcov[1][1]*gcov[2][3] - 
         gcov[0][1]*gcov[1][3]*gcov[2][0] - 
         gcov[0][3]*gcov[1][0]*gcov[2][1])/gdet;

    gcon[3][0] = gcon[0][3];
    gcon[3][1] = gcon[1][3];
    gcon[3][2] = gcon[2][3];

    gcon[3][3] =
        (gcov[0][0]*gcov[1][1]*gcov[2][2] + 
         gcov[0][1]*gcov[1][2]*gcov[2][0] + 
         gcov[0][2]*gcov[1][0]*gcov[2][1] - 
         gcov[0][0]*gcov[1][2]*gcov[2][1] - 
         gcov[0][1]*gcov[1][0]*gcov[2][2] - 
         gcov[0][2]*gcov[1][1]*gcov[2][0])/gdet;

}

void gDetCalc(REAL* gdet,
              const REAL gcov[NDIM][NDIM])
{
    *gdet = 
        gcov[0][0]*gcov[1][1]*gcov[2][2]*gcov[3][3] + 
        gcov[0][0]*gcov[1][2]*gcov[2][3]*gcov[3][1] + 
        gcov[0][0]*gcov[1][3]*gcov[2][1]*gcov[3][2] + 
        gcov[0][1]*gcov[1][0]*gcov[2][3]*gcov[3][2] + 
        gcov[0][1]*gcov[1][2]*gcov[2][0]*gcov[3][3] + 
        gcov[0][1]*gcov[1][3]*gcov[2][2]*gcov[3][0] + 
        gcov[0][2]*gcov[1][0]*gcov[2][1]*gcov[3][3] + 
        gcov[0][2]*gcov[1][1]*gcov[2][3]*gcov[3][0] + 
        gcov[0][2]*gcov[1][3]*gcov[2][0]*gcov[3][1] + 
        gcov[0][3]*gcov[1][0]*gcov[2][2]*gcov[3][1] + 
        gcov[0][3]*gcov[1][1]*gcov[2][0]*gcov[3][2] + 
        gcov[0][3]*gcov[1][2]*gcov[2][1]*gcov[3][0] - 
        gcov[0][0]*gcov[1][1]*gcov[2][3]*gcov[3][2] - 
        gcov[0][0]*gcov[1][2]*gcov[2][1]*gcov[3][3] - 
        gcov[0][0]*gcov[1][3]*gcov[2][2]*gcov[3][1] - 
        gcov[0][1]*gcov[1][0]*gcov[2][2]*gcov[3][3] - 
        gcov[0][1]*gcov[1][2]*gcov[2][3]*gcov[3][0] - 
        gcov[0][1]*gcov[1][3]*gcov[2][0]*gcov[3][2] - 
        gcov[0][2]*gcov[1][0]*gcov[2][3]*gcov[3][1] - 
        gcov[0][2]*gcov[1][1]*gcov[2][0]*gcov[3][3] - 
        gcov[0][2]*gcov[1][3]*gcov[2][1]*gcov[3][0] - 
        gcov[0][3]*gcov[1][0]*gcov[2][1]*gcov[3][2] - 
        gcov[0][3]*gcov[1][1]*gcov[2][2]*gcov[3][0] - 
        gcov[0][3]*gcov[1][2]*gcov[2][0]*gcov[3][1];
}

void alphaCalc(REAL* alpha,
               const REAL gcon[NDIM][NDIM])
{
    *alpha = 1./sqrt(-gcon[0][0]);
}

void gammaCalc(REAL* gamma,
               const REAL var[DOF],
               const REAL gcov[NDIM][NDIM])
{
    *gamma = 
        sqrt(1 + gcov[1][1]*var[U1]*var[U1] + 
                 gcov[2][2]*var[U2]*var[U2] + 
                 gcov[3][3]*var[U3]*var[U3] + 
              2*(gcov[1][2]*var[U1]*var[U2] + 
                 gcov[1][3]*var[U1]*var[U3] + 
                 gcov[2][3]*var[U2]*var[U3]));
}

void dgammaCalc_dt(REAL* dgamma_dt,
                   const REAL gamma,
                   const REAL var[DOF],
                   const REAL dvar_dt[DOF],
                   const REAL gcov[NDIM][NDIM])
{
    *dgamma_dt = 
        ((gcov[1][1]*var[U1]*dvar_dt[U1] + 
          gcov[2][2]*var[U2]*dvar_dt[U2] + 
          gcov[3][3]*var[U3]*dvar_dt[U3])+ 
         (gcov[1][2]*dvar_dt[U1]*var[U2] + 
          gcov[1][2]*var[U1]*dvar_dt[U2] + 
          gcov[1][3]*dvar_dt[U1]*var[U3] + 
          gcov[1][3]*var[U1]*dvar_dt[U3] + 
          gcov[2][3]*dvar_dt[U2]*var[U3] + 
          gcov[2][3]*var[U2]*dvar_dt[U3]))/gamma;
}

void uconCalc(REAL ucon[NDIM],
              const REAL gamma,
              const REAL alpha,
              const REAL var[DOF],
              const REAL gcon[NDIM][NDIM])
{
    ucon[0] = gamma/alpha;
    ucon[1] = var[U1] - gamma*gcon[0][1]*alpha;
    ucon[2] = var[U2] - gamma*gcon[0][2]*alpha;
    ucon[3] = var[U3] - gamma*gcon[0][3]*alpha;
}

void duconCalc_dt(REAL ducon_dt[NDIM],
                  const REAL dgamma_dt,
                  const REAL alpha,
                  const REAL dvar_dt[DOF],
                  const REAL gcon[NDIM][NDIM])
{
    ducon_dt[0] = dgamma_dt/alpha;
    ducon_dt[1] = dvar_dt[U1] - dgamma_dt*gcon[0][1]*alpha;
    ducon_dt[2] = dvar_dt[U2] - dgamma_dt*gcon[0][2]*alpha;
    ducon_dt[3] = dvar_dt[U3] - dgamma_dt*gcon[0][3]*alpha;
}

void covFromCon(REAL cov[NDIM],
                const REAL con[NDIM],
                const REAL gcov[NDIM][NDIM])
{
    cov[0] = gcov[0][0]*con[0] + gcov[0][1]*con[1] +
             gcov[0][2]*con[2] + gcov[0][3]*con[3];

    cov[1] = gcov[1][0]*con[0] + gcov[1][1]*con[1] +
             gcov[1][2]*con[2] + gcov[1][3]*con[3];

    cov[2] = gcov[2][0]*con[0] + gcov[2][1]*con[1] +
             gcov[2][2]*con[2] + gcov[2][3]*con[3];

    cov[3] = gcov[3][0]*con[0] + gcov[3][1]*con[1] +
             gcov[3][2]*con[2] + gcov[3][3]*con[3];
}

void conFromCov(REAL con[NDIM],
                const REAL cov[NDIM],
                const REAL gcon[NDIM][NDIM])
{
    con[0] = gcon[0][0]*cov[0] + gcon[0][1]*cov[1] +
             gcon[0][2]*cov[2] + gcon[0][3]*cov[3];

    con[1] = gcon[1][0]*cov[0] + gcon[1][1]*cov[1] +
             gcon[1][2]*cov[2] + gcon[1][3]*cov[3];

    con[2] = gcon[2][0]*cov[0] + gcon[2][1]*cov[1] +
             gcon[2][2]*cov[2] + gcon[2][3]*cov[3];

    con[3] = gcon[3][0]*cov[0] + gcon[3][1]*cov[1] +
             gcon[3][2]*cov[2] + gcon[3][3]*cov[3];

}

void conDotCov(REAL* ans,
               const REAL con[NDIM],
               const REAL cov[NDIM])
{
    *ans = con[0]*cov[0] + con[1]*cov[1] +
           con[2]*cov[2] + con[3]*cov[3];
}

void bconCalc(REAL bcon[NDIM],
              const REAL var[DOF],
              const REAL ucon[NDIM],
              const REAL ucov[NDIM])
{
    bcon[0] = var[B1]*ucov[1]+ var[B2]*ucov[2]+ var[B3]*ucov[3];
    
    bcon[1] = (var[B1] + bcon[0]*ucon[1])/ucon[0];
    bcon[2] = (var[B2] + bcon[0]*ucon[2])/ucon[0];
    bcon[3] = (var[B3] + bcon[0]*ucon[3])/ucon[0];
}

void dbconCalc_dt(REAL dbcon_dt[NDIM],
                  const REAL ucon[NDIM],
                  const REAL ducon_dt[NDIM],
                  const REAL ucov[NDIM],
                  const REAL ducov_dt[NDIM],
                  const REAL bcon[NDIM],
                  const REAL var[DOF],
                  const REAL dvar_dt[DOF])
{

    dbcon_dt[0] = 
        (dvar_dt[B1]*ucov[1] + dvar_dt[B2]*ucov[2] + dvar_dt[B3]*ucov[3] +
         var[B1]*ducov_dt[1] + var[B2]*ducov_dt[2] + var[B3]*ducov_dt[3]);

    dbcon_dt[1] =
        (-(var[B1] + bcon[0]*ucon[1])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B1] + bcon[0]*ducon_dt[1] + dbcon_dt[0]*ucon[1])/ucon[0]);

    dbcon_dt[2] =
        (-(var[B2] + bcon[0]*ucon[2])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B2] + bcon[0]*ducon_dt[2] + dbcon_dt[0]*ucon[2])/ucon[0]);

    dbcon_dt[3] =
        (-(var[B3] + bcon[0]*ucon[3])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B3] + bcon[0]*ducon_dt[3] + dbcon_dt[0]*ucon[3])/ucon[0]);
}

void bSqrCalc(REAL* bsqr,
              const REAL bcon[NDIM],
              const REAL bcov[NDIM])
{
    *bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] +
            bcon[2]*bcov[2] + bcon[3]*bcov[3];
}

void mhdCalc(REAL mhd[NDIM][NDIM],
             const REAL var[DOF],
             const REAL ucon[NDIM],
             const REAL ucov[NDIM],
             const REAL bcon[NDIM],
             const REAL bcov[NDIM])
{
    REAL P = (ADIABATIC_INDEX - 1.)*var[UU];
    REAL bsqr;
    bSqrCalc(&bsqr, bcon, bcov);
    
#define DELTA(mu, nu) (mu==nu ? 1 : 0)

    for (int mu=0; mu<NDIM; mu++)
        for (int nu=0; nu<NDIM; nu++) {
            mhd[mu][nu] = (var[RHO] + var[UU] + P + bsqr)*ucon[mu]*ucov[nu] +
                          (P + 0.5*bsqr)*DELTA(mu, nu) - bcon[mu]*bcov[nu];
        }

#undef DELTA
}

void addSources(REAL dU_dt[DOF],
                const REAL ucon[NDIM],
                const REAL ucov[NDIM],
                const REAL bcon[NDIM],
                const REAL bcov[NDIM],
                const REAL gcon[NDIM][NDIM],
                const REAL gcov[NDIM][NDIM],
                const REAL mhd[NDIM][NDIM],
                const REAL var[DOF],
                const REAL g,
                const REAL X1, const REAL X2)
{
    
    REAL gcovh[NDIM][NDIM], gcovl[NDIM][NDIM];
    REAL conntmp[NDIM][NDIM][NDIM], conn[NDIM][NDIM][NDIM];
    REAL Xl[NDIM], Xh[NDIM];

    for (int k = 0; k < NDIM; k++) {
        Xl[0] = 0.; Xl[1] = X1; Xl[2] = X2; Xl[3] = 0.;
        Xh[0] = 0.; Xh[1] = X1; Xh[2] = X2; Xh[3] = 0.;
        Xl[k] = Xl[k] - EPS;
        Xh[k] = Xh[k] + EPS;
        gCovCalc(gcovh, Xh[1], Xh[2]);
        gCovCalc(gcovl, Xl[1], Xl[2]);
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				conn[i][j][k] = (gcovh[i][j] - gcovl[i][j])/(Xh[k] - Xl[k]);
            }
        }
    }

	/* now rearrange to find \Gamma_{ijk} */
	for (int i = 0; i < NDIM; i++)
		for (int j = 0; j < NDIM; j++)
			for (int k = 0; k < NDIM; k++)
			    conntmp[i][j][k] =
				    0.5 * (conn[j][i][k] + conn[k][i][j] -
					   conn[k][j][i]);

	/* finally, raise index */
	for (int i = 0; i < NDIM; i++)
		for (int j = 0; j < NDIM; j++)
			for (int k = 0; k < NDIM; k++) {
				conn[i][j][k] = 0.;
				for (int l = 0; l < NDIM; l++)
					conn[i][j][k] += gcon[i][l]*conntmp[l][j][k];
			}
    
    for (int j=0; j<NDIM; j++)
        for (int k=0; k<NDIM; k++) {
            dU_dt[UU] = dU_dt[UU] - g*(mhd[j][k]*conn[k][0][j]);
            dU_dt[U1] = dU_dt[U1] - g*(mhd[j][k]*conn[k][1][j]);
            dU_dt[U2] = dU_dt[U2] - g*(mhd[j][k]*conn[k][2][j]);
            dU_dt[U3] = dU_dt[U3] - g*(mhd[j][k]*conn[k][3][j]);

        }
}

void ComputeFluxAndU(REAL flux[DOF],
                     REAL U[DOF],
                     const REAL ucon[NDIM],
                     const REAL ucov[NDIM],
                     const REAL bcon[NDIM],
                     const REAL bcov[NDIM],
                     const REAL gcon[NDIM][NDIM],
                     const REAL gcov[NDIM][NDIM],
                     const REAL mhd[NDIM][NDIM],
                     const REAL var[DOF],
                     const REAL g,
                     const int dir)
{
    flux[RHO] = g*var[RHO]*ucon[dir];

    flux[UU] = g*mhd[dir][0];
    flux[U1] = g*mhd[dir][1];
    flux[U2] = g*mhd[dir][2];
    flux[U3] = g*mhd[dir][3];

    flux[B1] = g*(bcon[1]*ucon[dir] - bcon[dir]*ucon[1]);
    flux[B2] = g*(bcon[2]*ucon[dir] - bcon[dir]*ucon[2]);
    flux[B3] = g*(bcon[3]*ucon[dir] - bcon[dir]*ucon[3]);

    U[RHO] = g*var[RHO]*ucon[0];

    U[UU] = g*mhd[0][0];
    U[U1] = g*mhd[0][1];
    U[U2] = g*mhd[0][2];
    U[U3] = g*mhd[0][3];

    U[B1] = g*(bcon[1]*ucon[0] - bcon[0]*ucon[1]);
    U[B2] = g*(bcon[2]*ucon[0] - bcon[0]*ucon[2]);
    U[B3] = g*(bcon[3]*ucon[0] - bcon[0]*ucon[3]);

}

void ComputedU_dt(REAL dU_dt[DOF],
                  const REAL ucon[NDIM],
                  const REAL ducon_dt[NDIM],
                  const REAL ucov[NDIM],
                  const REAL ducov_dt[NDIM],
                  const REAL bcon[NDIM],
                  const REAL dbcon_dt[NDIM],
                  const REAL bcov[NDIM],
                  const REAL dbcov_dt[NDIM],
                  const REAL gcon[NDIM][NDIM],
                  const REAL gcov[NDIM][NDIM],
                  const REAL var[DOF],
                  const REAL dvar_dt[DOF],
                  const REAL gamma,
                  const REAL dgamma_dt,
                  const REAL alpha,
                  const REAL g)
{
    REAL P, dP_dt, bsqr, dbsqr_dt, tmp1, dtmp1_dt, tmp2, dtmp2_dt;

    bSqrCalc(&bsqr, bcon, bcov);

    dbsqr_dt = bcon[0]*dbcov_dt[0] + dbcon_dt[0]*bcov[0] +
               bcon[1]*dbcov_dt[1] + dbcon_dt[1]*bcov[1] +
               bcon[2]*dbcov_dt[2] + dbcon_dt[2]*bcov[2] +
               bcon[3]*dbcov_dt[3] + dbcon_dt[3]*bcov[3];

    P = (ADIABATIC_INDEX-1.)*var[UU];
    dP_dt = (ADIABATIC_INDEX-1.)*dvar_dt[UU];

    tmp1 = P + var[RHO] + var[UU] + bsqr;
    dtmp1_dt = dP_dt + dvar_dt[RHO] + dvar_dt[UU] + dbsqr_dt;

    tmp2 = P + 0.5*bsqr;
    dtmp2_dt = dP_dt + 0.5*dbsqr_dt;

    dU_dt[RHO] = g*(dvar_dt[RHO]*ucon[0] + var[RHO]*ducon_dt[0]);

    dU_dt[UU] = g*(dtmp1_dt*ucon[0]*ucov[0] +
                   tmp1*(ducon_dt[0]*ucov[0] + ucon[0]*ducov_dt[0]) +
                   dtmp2_dt - dbcon_dt[0]*bcov[0] - bcon[0]*dbcov_dt[0]);

    dU_dt[U1] = g*(dtmp1_dt*ucon[0]*ucov[1] +
                   tmp1*(ducon_dt[0]*ucov[1] + ucon[0]*ducov_dt[1]) -
                   dbcon_dt[0]*bcov[1] - bcon[0]*dbcov_dt[1]);

    dU_dt[U2] = g*(dtmp1_dt*ucon[0]*ucov[2] +
                   tmp1*(ducon_dt[0]*ucov[2] + ucon[0]*ducov_dt[2]) -
                   dbcon_dt[0]*bcov[2] - bcon[0]*dbcov_dt[2]);

    dU_dt[U3] = g*(dtmp1_dt*ucon[0]*ucov[3] +
                   tmp1*(ducon_dt[0]*ucov[3] + ucon[0]*ducov_dt[3]) -
                   dbcon_dt[0]*bcov[3] - bcon[0]*dbcov_dt[3]);

    dU_dt[B1] = g*dvar_dt[B1];
    dU_dt[B2] = g*dvar_dt[B2];
    dU_dt[B3] = g*dvar_dt[B3];

}
                  
void VChar(REAL* vmin, REAL* vmax,
           const REAL ucon[DOF], const REAL ucov[DOF],
           const REAL bsqr,
           const REAL gcon[NDIM][NDIM],
           const REAL var[DOF], const int dir)
{
    REAL Acov[NDIM], Acon[NDIM], Bcov[NDIM], Bcon[NDIM];
    REAL Asqr, Bsqr, vasqr, cssqr, cmsqr, Adotu, Bdotu, AdotB;
    REAL A, B, C, discr;
    REAL vp, vm;

    vasqr = bsqr/(bsqr + var[RHO] + ADIABATIC_INDEX*var[UU]);
    cssqr = (ADIABATIC_INDEX)*(ADIABATIC_INDEX-1)*var[UU]/(var[RHO] +
             ADIABATIC_INDEX*var[UU]);

    cmsqr = cssqr + vasqr - cssqr*vasqr;
    
    for (int mu=0; mu<NDIM; mu++)
        Acov[mu] = 0.;
    Acov[dir] = 1.;
    conFromCov(Acon, Acov, gcon);

    for (int mu=0; mu<NDIM; mu++)
        Bcov[mu] = 0.;
    Bcov[0] = 1.;
    conFromCov(Bcon, Bcov, gcon);

    conDotCov(&Asqr, Acon, Acov);
    conDotCov(&Bsqr, Bcon, Bcov);
    conDotCov(&Adotu, ucon, Acov);
    conDotCov(&Bdotu, ucon, Bcov);
    conDotCov(&AdotB, Acon, Bcov);

    A = (Bdotu*Bdotu) - (Bsqr + Bdotu*Bdotu)*cmsqr;
    B = 2.*(Adotu*Bdotu - (AdotB + Adotu*Bdotu)*cmsqr);
    C = Adotu*Adotu - (Asqr + Adotu*Adotu)*cmsqr;

    discr = sqrt(B*B - 4.*A*C);

    vp = -(-B + discr)/(2.*A);
    vm = -(-B - discr)/(2.*A);

    if (vp>vm) {
        *vmax = vp;
        *vmin = vm;
    } else {
        *vmax = vm;
        *vmin = vp;
    }
}

void setFloor(REAL var[DOF], 
              const REAL gamma,
              REAL X1, REAL X2)
{
    REAL r, theta;
    BLCoords(&r, &theta, X1, X2);

    REAL rhoFloor = RHO_MIN*pow(r, -1.5);
    REAL uFloor = U_MIN*pow(r, -2.5);

    if (rhoFloor < RHO_MIN_LIMIT)
        rhoFloor = RHO_MIN_LIMIT;

    if (uFloor < U_MIN_LIMIT)
        uFloor = U_MIN_LIMIT;
    
    if (var[RHO] < rhoFloor)
        var[RHO] = rhoFloor;

    if (var[UU] < uFloor)
        var[UU] = uFloor;

    if (gamma > GAMMA_MAX) {
        REAL f = sqrt((GAMMA_MAX*GAMMA_MAX - 1.) /
				      (gamma * gamma - 1.));

        var[U1] = f*var[U1];
        var[U2] = f*var[U2];
        var[U3] = f*var[U3];

    }
}

//#define LEFT -1
//#define RIGHT 1
//
//REAL SlopeLim(REAL y1, REAL y2, REAL y3)
//{
//    REAL Dqm = 2. * (y2 - y1);
//	REAL Dqp = 2. * (y3 - y2);
//	REAL Dqc = 0.5 * (y3 - y1);
//	REAL s = Dqm * Dqp;
//	if (s <= 0.) {
//		return 0.;
//    }
//	else {
//		if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
//			return (Dqm);
//		else if (fabs(Dqp) < fabs(Dqc))
//			return (Dqp);
//		else
//			return (Dqc);
//	}
//
//}
//
//void ReconstructX1(const REAL* primTile, 
//                   const int iTile, const int jTile,
//                   REAL primEdge[DOF],
//                   const int dir)
//{
//    REAL left, mid, right, slope;
//    for (int var=0; var<DOF; var++) {
//        left = primTile[INDEX_LOCAL(iTile-1, jTile, var)];
//        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
//        right = primTile[INDEX_LOCAL(iTile+1, jTile, var)];
//        
//        slope = SlopeLim(left, mid, right);
//
//        primEdge[var] = mid + dir*0.5*slope;
//    }
//}
//
//void ReconstructX2(const REAL* primTile, 
//                   const int iTile, const int jTile,
//                   REAL primEdge[DOF],
//                   const int dir)
//{
//    REAL left, mid, right, slope;
//    for (int var=0; var<DOF; var++) {
//        left = primTile[INDEX_LOCAL(iTile, jTile-1, var)];
//        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
//        right = primTile[INDEX_LOCAL(iTile, jTile+1, var)];
//        
//        slope = SlopeLim(left, mid, right);
//
//        primEdge[var] = mid + dir*0.5*slope;
//    }
//
//}
//
//void RiemannSolver(const REAL fluxL[DOF],
//                   const REAL fluxR[DOF],
//                   const REAL uL[DOF],
//                   const REAL uR[DOF],
//                   REAL flux[DOF])
//{
//    for (int var=0; var<DOF; var++) {
//        flux[var] = 0.5*(fluxL[var] + fluxR[var] +
//                         (uL[var] - uR[var]));
//    }
//}
////Boyer Linquist coordinates
//#define R_BL(X1,X2) (exp(X1) + R0)
//#define THETA_BL(X1,X2) (M_PI*(X2) + ((1 - H_SLOPE)/2.)*sin(2.*M_PI*(X2)))
//
////Kerr Schild metric
////Parameters in the Kerr Schild metric
//#define RHO_SQR(X1,X2) \
//    (R_BL(X1,X2)*R_BL(X1,X2) + \
//     A_SPIN*A_SPIN*cos(THETA_BL(X1,X2))*cos(THETA_BL(X1,X2)))
//
//#define R_FACTOR(X1,X2) (R_BL(X1,X2) - R0)
//#define H_FACTOR(X1,X2) (M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*(X2)))
//
//#define GCOV00(X1,X2) (-1. + 2.*R_BL(X1,X2)/RHO_SQR(X1,X2))
//#define GCOV01(X1,X2) ((2.*R_BL(X1,X2)/RHO_SQR(X1,X2))*R_FACTOR(X1,X2))
//#define GCOV02(X1,X2) (0.)
//#define GCOV03(X1,X2) \
//    (-2.*A_SPIN*R_BL(X1,X2)*sin(THETA_BL(X1,X2))*sin(THETA_BL(X1,X2))\
//     /RHO_SQR(X1,X2) )
//
//#define GCOV10(X1,X2) (GCOV01(X1,X2))
//#define GCOV11(X1,X2) \
//    ( (1. + 2.*R_BL(X1,X2)/RHO_SQR(X1,X2))*R_FACTOR(X1,X2)*R_FACTOR(X1,X2) )
//
//#define GCOV12(X1,X2) (0.)
//#define GCOV13(X1,X2) \
//    (-A_SPIN*sin(THETA_BL(X1,X2))*sin(THETA_BL(X1,X2))* \
//     (1.+2.*R_BL(X1,X2)/RHO_SQR(X1,X2))*R_FACTOR(X1,X2))
//
//#define GCOV20(X1,X2) (GCOV02(X1,X2))
//#define GCOV21(X1,X2) (GCOV12(X1,X2))
//#define GCOV22(X1,X2) (RHO_SQR(X1,X2)*H_FACTOR(X1,X2)*H_FACTOR(X1,X2))
//#define GCOV23(X1,X2) (0.)
//#define GCOV30(X1,X2) (GCOV03(X1,X2))
//#define GCOV31(X1,X2) (GCOV13(X1,X2))
//#define GCOV32(X1,X2) (GCOV23(X1,X2))
//#define GCOV33(X1,X2) \
//    (sin(THETA_BL(X1,X2))*sin(THETA_BL(X1,X2))*(RHO_SQR(X1,X2) + \
//     A_SPIN*A_SPIN*sin(THETA_BL(X1,X2))*sin(THETA_BL(X1,X2))*\
//     (1. + 2.*R_BL(X1,X2)/RHO_SQR(X1,X2)) ) )
//
//
//#define GCOV00(X1,X2) (-1.)
//#define GCOV01(X1,X2) (0.)
//#define GCOV02(X1,X2) (0.)
//#define GCOV03(X1,X2) (0.)
//
//#define GCOV10(X1,X2) (0.)
//#define GCOV11(X1,X2) (1.)
//#define GCOV12(X1,X2) (0.)
//#define GCOV13(X1,X2) (0.)
//
//#define GCOV20(X1,X2) (0.)
//#define GCOV21(X1,X2) (0.)
//#define GCOV22(X1,X2) (1.)
//#define GCOV23(X1,X2) (0.)
//
//#define GCOV30(X1,X2) (0.)
//#define GCOV31(X1,X2) (0.)
//#define GCOV32(X1,X2) (0.)
//#define GCOV33(X1,X2) (1.)

//#define GDET(X1,X2) \
//    (GCOV00(X1,X2)*GCOV11(X1,X2)*GCOV22(X1,X2)*GCOV33(X1,X2) + \
//     GCOV00(X1,X2)*GCOV12(X1,X2)*GCOV23(X1,X2)*GCOV31(X1,X2) + \
//     GCOV00(X1,X2)*GCOV13(X1,X2)*GCOV21(X1,X2)*GCOV32(X1,X2) + \
//     GCOV01(X1,X2)*GCOV10(X1,X2)*GCOV23(X1,X2)*GCOV32(X1,X2) + \
//     GCOV01(X1,X2)*GCOV12(X1,X2)*GCOV20(X1,X2)*GCOV33(X1,X2) + \
//     GCOV01(X1,X2)*GCOV13(X1,X2)*GCOV22(X1,X2)*GCOV30(X1,X2) + \
//     GCOV02(X1,X2)*GCOV10(X1,X2)*GCOV21(X1,X2)*GCOV33(X1,X2) + \
//     GCOV02(X1,X2)*GCOV11(X1,X2)*GCOV23(X1,X2)*GCOV30(X1,X2) + \
//     GCOV02(X1,X2)*GCOV13(X1,X2)*GCOV20(X1,X2)*GCOV31(X1,X2) + \
//     GCOV03(X1,X2)*GCOV10(X1,X2)*GCOV22(X1,X2)*GCOV31(X1,X2) + \
//     GCOV03(X1,X2)*GCOV11(X1,X2)*GCOV20(X1,X2)*GCOV32(X1,X2) + \
//     GCOV03(X1,X2)*GCOV12(X1,X2)*GCOV21(X1,X2)*GCOV30(X1,X2) - \
//     GCOV00(X1,X2)*GCOV11(X1,X2)*GCOV23(X1,X2)*GCOV32(X1,X2) - \
//     GCOV00(X1,X2)*GCOV12(X1,X2)*GCOV21(X1,X2)*GCOV33(X1,X2) - \
//     GCOV00(X1,X2)*GCOV13(X1,X2)*GCOV22(X1,X2)*GCOV31(X1,X2) - \
//     GCOV01(X1,X2)*GCOV10(X1,X2)*GCOV22(X1,X2)*GCOV33(X1,X2) - \
//     GCOV01(X1,X2)*GCOV12(X1,X2)*GCOV23(X1,X2)*GCOV30(X1,X2) - \
//     GCOV01(X1,X2)*GCOV13(X1,X2)*GCOV20(X1,X2)*GCOV32(X1,X2) - \
//     GCOV02(X1,X2)*GCOV10(X1,X2)*GCOV23(X1,X2)*GCOV31(X1,X2) - \
//     GCOV02(X1,X2)*GCOV11(X1,X2)*GCOV20(X1,X2)*GCOV33(X1,X2) - \
//     GCOV02(X1,X2)*GCOV13(X1,X2)*GCOV21(X1,X2)*GCOV30(X1,X2) - \
//     GCOV03(X1,X2)*GCOV10(X1,X2)*GCOV21(X1,X2)*GCOV32(X1,X2) - \
//     GCOV03(X1,X2)*GCOV11(X1,X2)*GCOV22(X1,X2)*GCOV30(X1,X2) - \
//     GCOV03(X1,X2)*GCOV12(X1,X2)*GCOV20(X1,X2)*GCOV31(X1,X2))
//
//#define GCON00(X1,X2,gDet) \
//    ((GCOV11(X1,X2)*GCOV22(X1,X2)*GCOV33(X1,X2) + \
//      GCOV12(X1,X2)*GCOV23(X1,X2)*GCOV31(X1,X2) + \
//      GCOV13(X1,X2)*GCOV21(X1,X2)*GCOV32(X1,X2) - \
//      GCOV11(X1,X2)*GCOV23(X1,X2)*GCOV32(X1,X2) - \
//      GCOV12(X1,X2)*GCOV21(X1,X2)*GCOV33(X1,X2) - \
//      GCOV13(X1,X2)*GCOV22(X1,X2)*GCOV31(X1,X2))/gDet) 
//
//#define GCON01(X1,X2,gDet) \
//    ((GCOV01(X1,X2)*GCOV23(X1,X2)*GCOV32(X1,X2) + \
//      GCOV02(X1,X2)*GCOV21(X1,X2)*GCOV33(X1,X2) + \
//      GCOV03(X1,X2)*GCOV22(X1,X2)*GCOV31(X1,X2) - \
//      GCOV01(X1,X2)*GCOV22(X1,X2)*GCOV33(X1,X2) - \
//      GCOV02(X1,X2)*GCOV23(X1,X2)*GCOV31(X1,X2) - \
//      GCOV03(X1,X2)*GCOV21(X1,X2)*GCOV32(X1,X2))/gDet) 
//
//#define GCON02(X1,X2,gDet) \
//    ((GCOV01(X1,X2)*GCOV12(X1,X2)*GCOV33(X1,X2) + \
//      GCOV02(X1,X2)*GCOV13(X1,X2)*GCOV31(X1,X2) + \
//      GCOV03(X1,X2)*GCOV11(X1,X2)*GCOV32(X1,X2) - \
//      GCOV01(X1,X2)*GCOV13(X1,X2)*GCOV32(X1,X2) - \
//      GCOV02(X1,X2)*GCOV11(X1,X2)*GCOV33(X1,X2) - \
//      GCOV03(X1,X2)*GCOV12(X1,X2)*GCOV31(X1,X2))/gDet) 
//
//#define GCON03(X1,X2,gDet) \
//    ((GCOV01(X1,X2)*GCOV13(X1,X2)*GCOV22(X1,X2) + \
//      GCOV02(X1,X2)*GCOV11(X1,X2)*GCOV23(X1,X2) + \
//      GCOV03(X1,X2)*GCOV12(X1,X2)*GCOV21(X1,X2) - \
//      GCOV01(X1,X2)*GCOV12(X1,X2)*GCOV23(X1,X2) - \
//      GCOV02(X1,X2)*GCOV13(X1,X2)*GCOV21(X1,X2) - \
//      GCOV03(X1,X2)*GCOV11(X1,X2)*GCOV22(X1,X2))/gDet) 
//
//#define GCON10(X1,X2,gDet) (GCON01(X1,X2,gDet))
//    
//#define GCON11(X1,X2,gDet) \
//    ((GCOV00(X1,X2)*GCOV22(X1,X2)*GCOV33(X1,X2) + \
//      GCOV02(X1,X2)*GCOV23(X1,X2)*GCOV30(X1,X2) + \
//      GCOV03(X1,X2)*GCOV20(X1,X2)*GCOV32(X1,X2) - \
//      GCOV00(X1,X2)*GCOV23(X1,X2)*GCOV32(X1,X2) - \
//      GCOV02(X1,X2)*GCOV20(X1,X2)*GCOV33(X1,X2) - \
//      GCOV03(X1,X2)*GCOV22(X1,X2)*GCOV30(X1,X2))/gDet) 
//
//#define GCON12(X1,X2,gDet) \
//    ((GCOV00(X1,X2)*GCOV13(X1,X2)*GCOV32(X1,X2) + \
//      GCOV02(X1,X2)*GCOV10(X1,X2)*GCOV33(X1,X2) + \
//      GCOV03(X1,X2)*GCOV12(X1,X2)*GCOV30(X1,X2) - \
//      GCOV00(X1,X2)*GCOV12(X1,X2)*GCOV33(X1,X2) - \
//      GCOV02(X1,X2)*GCOV13(X1,X2)*GCOV30(X1,X2) - \
//      GCOV03(X1,X2)*GCOV10(X1,X2)*GCOV32(X1,X2))/gDet) 
//
//#define GCON13(X1,X2,gDet) \
//    ((GCOV00(X1,X2)*GCOV12(X1,X2)*GCOV23(X1,X2) + \
//      GCOV02(X1,X2)*GCOV13(X1,X2)*GCOV20(X1,X2) + \
//      GCOV03(X1,X2)*GCOV10(X1,X2)*GCOV22(X1,X2) - \
//      GCOV00(X1,X2)*GCOV13(X1,X2)*GCOV22(X1,X2) - \
//      GCOV02(X1,X2)*GCOV10(X1,X2)*GCOV23(X1,X2) - \
//      GCOV03(X1,X2)*GCOV12(X1,X2)*GCOV20(X1,X2))/gDet) 
//
//#define GCON20(X1,X2,gDet) (GCON02(X1,X2,gDet))
//#define GCON21(X1,X2,gDet) (GCON12(X1,X2,gDet))
//
//#define GCON22(X1,X2,gDet) \
//    ((GCOV00(X1,X2)*GCOV11(X1,X2)*GCOV33(X1,X2) + \
//      GCOV01(X1,X2)*GCOV13(X1,X2)*GCOV30(X1,X2) + \
//      GCOV03(X1,X2)*GCOV10(X1,X2)*GCOV31(X1,X2) - \
//      GCOV00(X1,X2)*GCOV13(X1,X2)*GCOV31(X1,X2) - \
//      GCOV01(X1,X2)*GCOV10(X1,X2)*GCOV33(X1,X2) - \
//      GCOV03(X1,X2)*GCOV11(X1,X2)*GCOV30(X1,X2))/gDet)
//
//#define GCON23(X1,X2,gDet) \
//    ((GCOV00(X1,X2)*GCOV13(X1,X2)*GCOV21(X1,X2) + \
//      GCOV01(X1,X2)*GCOV10(X1,X2)*GCOV23(X1,X2) + \
//      GCOV03(X1,X2)*GCOV11(X1,X2)*GCOV20(X1,X2) - \
//      GCOV00(X1,X2)*GCOV11(X1,X2)*GCOV23(X1,X2) - \
//      GCOV01(X1,X2)*GCOV13(X1,X2)*GCOV20(X1,X2) - \
//      GCOV03(X1,X2)*GCOV10(X1,X2)*GCOV21(X1,X2))/gDet)
//
//#define GCON30(X1,X2,gDet) (GCON03(X1,X2,gDet))
//#define GCON31(X1,X2,gDet) (GCON13(X1,X2,gDet))
//#define GCON32(X1,X2,gDet) (GCON23(X1,X2,gDet))
//
//#define GCON33(X1,X2,gDet) \
//    ((GCOV00(X1,X2)*GCOV11(X1,X2)*GCOV22(X1,X2) + \
//      GCOV01(X1,X2)*GCOV12(X1,X2)*GCOV20(X1,X2) + \
//      GCOV02(X1,X2)*GCOV10(X1,X2)*GCOV21(X1,X2) - \
//      GCOV00(X1,X2)*GCOV12(X1,X2)*GCOV21(X1,X2) - \
//      GCOV01(X1,X2)*GCOV10(X1,X2)*GCOV22(X1,X2) - \
//      GCOV02(X1,X2)*GCOV11(X1,X2)*GCOV20(X1,X2))/gDet)
//
//#define ALPHA (1./sqrt(-GCON00(X1,X2,gDet)))
//
//#define GAMMA (sqrt(1 + GCOV11(X1,X2)*var[U1]*var[U1] + \
//                        GCOV22(X1,X2)*var[U2]*var[U2] + \
//                        GCOV33(X1,X2)*var[U3]*var[U3] + \
//                     2*(GCOV12(X1,X2)*var[U1]*var[U2] + \
//                        GCOV13(X1,X2)*var[U1]*var[U3] + \
//                        GCOV23(X1,X2)*var[U2]*var[U3])) )
//
//#define DGAMMA_DT ( ((GCOV11(X1,X2)*var[U1]*dvar_dt[U1] + \
//                      GCOV22(X1,X2)*var[U2]*dvar_dt[U2] + \
//                      GCOV33(X1,X2)*var[U3]*dvar_dt[U3]) + \
//                     (GCOV12(X1,X2)*dvar_dt[U1]*var[U2] + \
//                      GCOV12(X1,X2)*var[U1]*dvar_dt[U2] + \
//                      GCOV13(X1,X2)*dvar_dt[U1]*var[U3] + \
//                      GCOV13(X1,X2)*var[U1]*dvar_dt[U3] + \
//                      GCOV23(X1,X2)*dvar_dt[U2]*var[U3] + \
//                      GCOV23(X1,X2)*var[U2]*dvar_dt[U3]))/GAMMA)
//
//#define UCON0 (gamma/alpha)
//#define UCON1 (var[U1] - gamma*GCON01(X1,X2,gDet)*alpha)
//#define UCON2 (var[U2] - gamma*GCON02(X1,X2,gDet)*alpha)
//#define UCON3 (var[U3] - gamma*GCON03(X1,X2,gDet)*alpha)
//
//#define DUCON0_DT (dgamma_dt/alpha)
//#define DUCON1_DT (dvar_dt[U1] - dgamma_dt*GCON01(X1,X2)*alpha)
//#define DUCON2_DT (dvar_dt[U2] - dgamma_dt*GCON02(X1,X2)*alpha)
//#define DUCON3_DT (dvar_dt[U3] - dgamma_dt*GCON03(X1,X2)*alpha)
//
//#define UCOV0 (GCOV00(X1,X2)*UCON0 + GCOV01(X1,X2)*UCON1 + \
//               GCOV02(X1,X2)*UCON2 + GCOV03(X1,X2)*UCON3 )
//#define UCOV1 (GCOV10(X1,X2)*UCON0 + GCOV11(X1,X2)*UCON1 + \
//               GCOV12(X1,X2)*UCON2 + GCOV13(X1,X2)*UCON3 )
//#define UCOV2 (GCOV20(X1,X2)*UCON0 + GCOV21(X1,X2)*UCON1 + \
//               GCOV22(X1,X2)*UCON2 + GCOV23(X1,X2)*UCON3 )
//#define UCOV3 (GCOV30(X1,X2)*UCON0 + GCOV31(X1,X2)*UCON1 + \
//               GCOV32(X1,X2)*UCON2 + GCOV33(X1,X2)*UCON3 )
//
//#define DUCOV0_DT (GCOV00(X1,X2)*DUCON0_DT + GCOV01(X1,X2)*DUCON1_DT + \
//                   GCOV02(X1,X2)*DUCON2_DT + GCOV03(X1,X2)*DUCON3_DT )
//#define DUCOV1_DT (GCOV10(X1,X2)*DUCON0_DT + GCOV11(X1,X2)*DUCON1_DT + \
//                   GCOV12(X1,X2)*DUCON2_DT + GCOV13(X1,X2)*DUCON3_DT )
//#define DUCOV2_DT (GCOV20(X1,X2)*DUCON0_DT + GCOV21(X1,X2)*DUCON1_DT + \
//                   GCOV22(X1,X2)*DUCON2_DT + GCOV23(X1,X2)*DUCON3_DT )
//#define DUCOV3_DT (GCOV30(X1,X2)*DUCON0_DT + GCOV31(X1,X2)*DUCON1_DT + \
//                   GCOV32(X1,X2)*DUCON2_DT + GCOV33(X1,X2)*DUCON3_DT )
//
//#define BCON0 (var[B1]*UCOV1+ var[B2]*UCOV2+ var[B3]*UCOV3)
//
//#define BCON1 ((var[B1] + BCON0*UCON1)/UCON0)
//#define BCON2 ((var[B2] + BCON0*UCON2)/UCON0)
//#define BCON3 ((var[B3] + BCON0*UCON3)/UCON0)
//
//#define DBCON0_DT (dvar_dt[B1]*UCOV1 + dvar_dt[B2]*UCOV2 + dvar_dt[B3]*UCOV3 +\
//                   var[B1]*DUCOV1_DT + var[B2]*DUCOV2_DT + var[B3]*DUCOV3_DT)
//
//#define DBCON1_DT (-(var[B1] + BCON0*UCON1)*DUCON0_DT/(UCON0*UCON0) + \
//                   (dvar_dt[B1] + BCON0*DUCON1_DT + DBCON0_DT*UCON1)/UCON0)
//
//#define DBCON2_DT (-(var[B2] + BCON0*UCON2)*DUCON0_DT/(UCON0*UCON0) + \
//                   (dvar_dt[B2] + BCON0*DUCON2_DT + DBCON0_DT*UCON2)/UCON0)
//
//#define DBCON3_DT (-(var[B3] + BCON0*UCON3)*DUCON0_DT/(UCON0*UCON0) + \
//                   (dvar_dt[B3] + BCON0*DUCON3_DT + DBCON0_DT*UCON3)/UCON0)
//
//#define BCOV0 (GCOV00(X1,X2)*BCON0 + GCOV01(X1,X2)*BCON1 + \
//               GCOV02(X1,X2)*BCON2 + GCOV03(X1,X2)*BCON3 )
//#define BCOV1 (GCOV10(X1,X2)*BCON0 + GCOV11(X1,X2)*BCON1 + \
//               GCOV12(X1,X2)*BCON2 + GCOV13(X1,X2)*BCON3 )
//#define BCOV2 (GCOV20(X1,X2)*BCON0 + GCOV21(X1,X2)*BCON1 + \
//               GCOV22(X1,X2)*BCON2 + GCOV23(X1,X2)*BCON3 )
//#define BCOV3 (GCOV30(X1,X2)*BCON0 + GCOV31(X1,X2)*BCON1 + \
//               GCOV32(X1,X2)*BCON2 + GCOV33(X1,X2)*BCON3 )
//
//#define DBCOV0_DT (GCOV00(X1,X2)*DBCON0_DT + GCOV01(X1,X2)*DBCON1_DT + \
//                   GCOV02(X1,X2)*DBCON2_DT + GCOV03(X1,X2)*DBCON3_DT )
//#define DBCOV1_DT (GCOV10(X1,X2)*DBCON0_DT + GCOV11(X1,X2)*DBCON1_DT + \
//                   GCOV12(X1,X2)*DBCON2_DT + GCOV13(X1,X2)*DBCON3_DT )
//#define DBCOV2_DT (GCOV20(X1,X2)*DBCON0_DT + GCOV21(X1,X2)*DBCON1_DT + \
//                   GCOV22(X1,X2)*DBCON2_DT + GCOV23(X1,X2)*DBCON3_DT )
//#define DBCOV3_DT (GCOV30(X1,X2)*DBCON0_DT + GCOV31(X1,X2)*DBCON1_DT + \
//                   GCOV32(X1,X2)*DBCON2_DT + GCOV33(X1,X2)*DBCON3_DT )
//
//#define BSQR (BCON0*BCOV0 + BCON1*BCOV1 + \
//              BCON2*BCOV2 + BCON3*BCOV3 )
//
//#define DBSQR_DT (BCON0*DBCOV0_DT + DBCON0_DT*BCOV0 + \
//                  BCON1*DBCOV1_DT + DBCON1_DT*BCOV1 + \
//                  BCON2*DBCOV2_DT + DBCON2_DT*BCOV2 + \
//                  BCON3*DBCOV3_DT + DBCON3_DT*BCOV3 )
//
//#define P ((ADIABATIC_INDEX - 1.)*var[UU])
//#define DP_DT ((ADIABATIC_INDEX - 1.)*dvar_dt[UU])
//
//#define MHD0_0 (tmp1*UCON0*UCOV0 + tmp2 - BCON0*BCOV0)
//#define MHD0_1 (tmp1*UCON0*UCOV1 - BCON0*BCOV1)
//#define MHD0_2 (tmp1*UCON0*UCOV2 - BCON0*BCOV2)
//#define MHD0_3 (tmp1*UCON0*UCOV3 - BCON0*BCOV3)
//
//#define MHD1_0 (tmp1*UCON1*UCOV0 - BCON1*BCOV0)
//#define MHD1_1 (tmp1*UCON1*UCOV1 + tmp2 - BCON1*BCOV1)
//#define MHD1_2 (tmp1*UCON1*UCOV2 - BCON1*BCOV2)
//#define MHD1_3 (tmp1*UCON1*UCOV3 - BCON1*BCOV3)
//
//#define MHD2_0 (tmp1*UCON2*UCOV0 - BCON2*BCOV0)
//#define MHD2_1 (tmp1*UCON2*UCOV1 - BCON2*BCOV1)
//#define MHD2_2 (tmp1*UCON2*UCOV2 + tmp2 - BCON2*BCOV2)
//#define MHD2_3 (tmp1*UCON2*UCOV3 - BCON2*BCOV3)
// 
//#define MHD3_0 (tmp1*UCON3*UCOV0 - BCON3*BCOV0)
//#define MHD3_1 (tmp1*UCON3*UCOV1 - BCON3*BCOV1)
//#define MHD3_2 (tmp1*UCON3*UCOV2 - BCON3*BCOV2)
//#define MHD3_3 (tmp1*UCON3*UCOV3 + tmp2 - BCON3*BCOV3)
