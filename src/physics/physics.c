#include "physics.h"

void setGamma(const struct geometry geom[ARRAY_ARGS 1],
              struct fluidElement elem[ARRAY_ARGS 1])
{
  elem->gamma = 
    sqrt(1 + geom->gCov[1][1]*elem->primVars[U1]*elem->primVars[U1]
           + geom->gCov[2][2]*elem->primVars[U2]*elem->primVars[U2] 
           + geom->gCov[3][3]*elem->primVars[U3]*elem->primVars[U3]

         + 2*(  geom->gCov[1][2]*elem->primVars[U1]*elem->primVars[U2]
              + geom->gCov[1][3]*elem->primVars[U1]*elem->primVars[U3] 
              + geom->gCov[2][3]*elem->primVars[U2]*elem->primVars[U3]
             )
        );
}

void setUCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1])
{
  elem->uCon[0] = elem->gamma/geom->alpha;

  for (int i=1; i<NDIM; i++)
  {
    elem->uCon[i] =  elem->primVars[UU+i] 
                   - elem->gamma*geom->gCon[0][i]*geom->alpha;
  }
}

void setBCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1])
{
  REAL uCov[NDIM];

  conToCov(elem->uCon, geom, uCov);

  elem->bCon[0] =   elem->primVars[B1]*uCov[1]
                  + elem->primVars[B2]*uCov[2] 
                  + elem->primVars[B3]*uCov[3];
  
  for (int i=1; i<NDIM; i++)
  {
    elem->bCon[i] = (  elem->primVars[U3+i] 
                     + elem->bCon[0]*elem->uCon[i]
                    )/elem->uCon[0];
  }
}

void setFluidElement(const REAL primVars[ARRAY_ARGS DOF],
                     const struct geometry geom[ARRAY_ARGS 1],
                     struct fluidElement elem[ARRAY_ARGS 1])
{
  for (int var=0; var<DOF; var++)
  {
    elem->primVars[var] = primVars[var];
  }

  /* Need to be set in exactly the following order because of the dependecy
     structure */
  setGamma(geom, elem);
  setUCon(geom, elem);
  setBCon(geom, elem);
  computeMoments(geom, elem);
}

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1])
{
  REAL pressure = (ADIABATIC_INDEX - 1.)*elem->primVars[UU];
  REAL bCov[NDIM], bSqr;
  
  conToCov(elem->bCon, geom, bCov);
  bSqr = covDotCon(bCov, elem->bCon);

  for (int mu=0; mu<NDIM; mu++)
  {
    elem->moments[N_UP(mu)] = elem->primVars[RHO]*elem->uCon[mu];

    for (int nu=0; nu<NDIM; nu++)
    {
      elem->moments[T_UP_UP(mu,nu)] =   
                          (  elem->primVars[RHO] + elem->primVars[UU]
                           + pressure + bSqr
                          )*elem->uCon[mu]*elem->uCon[nu]

                        + (pressure + 0.5*bSqr)*geom->gCon[mu][nu]

                        - elem->bCon[mu]*elem->bCon[nu];

    }
  }

}

void computeFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                   const struct geometry geom[ARRAY_ARGS 1],
                   const int dir,
                   REAL fluxes[ARRAY_ARGS DOF])
{
  REAL g = sqrt(-geom->gDet);

  fluxes[RHO] = g*elem->moments[N_UP(dir)];

  fluxes[UU] = g*elem->moments[T_UP_UP(dir, 0)];
  fluxes[U1] = g*elem->moments[T_UP_UP(dir, 1)];
  fluxes[U2] = g*elem->moments[T_UP_UP(dir, 2)];
  fluxes[U3] = g*elem->moments[T_UP_UP(dir, 3)];

  fluxes[B1] = g*(elem->bCon[1]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[1]);
  fluxes[B2] = g*(elem->bCon[2]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[2]);
  fluxes[B3] = g*(elem->bCon[3]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[3]);
}

/*  Returns sqrt(-gDet) * T^kappa^lamda * Gamma_lamda_nu_kappa for each nu.
 *  (Eqn (4) of HARM paper)
 * 
*/

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL X[ARRAY_ARGS NDIM],
                        REAL sourceTerms[ARRAY_ARGS DOF])
{
  REAL g = sqrt(-geom->gDet);

  for (int var=0; var<DOF; var++)
  {
    sourceTerms[var] = 0.;
  }

//  for (int nu=0; nu<NDIM; nu++)
//  {
//    sourceTerms[UU+nu] = 0.;
//
//    for (int kappa=0; kappa<NDIM; kappa++)
//    {
//      for (int lamda=0; lamda<NDIM; lamda++)
//      {
//        for (int mu=0; mu<NDIM; mu++)
//        {
//          sourceTerms[RHO+nu] +=   
//               g * elem->moments[T_UP_UP(kappa, lamda)]
//                 * gammaDownDownDown(lamda, nu, kappa, X);
//        }
//      }
//    }
//  }


  REAL conntmp[NDIM][NDIM][NDIM], conn[NDIM][NDIM][NDIM];
  REAL Xl[NDIM], Xh[NDIM];
  struct geometry geomh, geoml;
  
  for (int k = 0; k < NDIM; k++) {
    Xl[0] = 0.; Xl[1] = X[1]; Xl[2] = X[2]; Xl[3] = 0.;
    Xh[0] = 0.; Xh[1] = X[1]; Xh[2] = X[2]; Xh[3] = 0.;
    Xl[k] = Xl[k] - EPS;
    Xh[k] = Xh[k] + EPS;
    
    setGeometry(Xh, &geomh);
    setGeometry(Xl, &geoml);
    for (int i = 0; i < NDIM; i++) {
      for (int j = 0; j < NDIM; j++) {
        conn[i][j][k] = (geomh.gCov[i][j] - geoml.gCov[i][j])/(Xh[k] - Xl[k]);
      }
    }
  }
  /* now rearrange to find \Gamma_{ijk} */
  for (int i = 0; i < NDIM; i++)
    for (int j = 0; j < NDIM; j++)
      for (int k = 0; k < NDIM; k++)
        conntmp[i][j][k] =
          0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

  for (int i = 0; i < NDIM; i++)
    for (int j = 0; j < NDIM; j++)
      for (int k = 0; k < NDIM; k++) {
        conn[i][j][k] = 0.;
        for (int l = 0; l < NDIM; l++)
          conn[i][j][k] += geom->gCon[i][l]*conntmp[l][j][k];
      }

   for (int j=0; j<NDIM; j++)
    for (int k=0; k<NDIM; k++) {
      sourceTerms[UU] =  g*(elem->moments[T_UP_UP(j,k)]*conntmp[k][0][j]);
      sourceTerms[U1] =  g*(elem->moments[T_UP_UP(j,k)]*conntmp[k][1][j]);
      sourceTerms[U2] =  g*(elem->moments[T_UP_UP(j,k)]*conntmp[k][2][j]);
      sourceTerms[U3] =  g*(elem->moments[T_UP_UP(j,k)]*conntmp[k][3][j]);
    }

}
