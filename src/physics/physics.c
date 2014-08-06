#include "physics.h"

void setGamma(const struct geometry* restrict geom,
              struct fluidElement* restrict elem)
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

void setUCon(const struct geometry* restrict geom,
             struct fluidElement* restrict elem)
{
  elem->uCon[0] = elem->gamma/geom->alpha;

  for (int i=1; i<NDIM; i++)
  {
    elem->uCon[i] =   elem->primVars[UU+i] 
                   - elem->gamma*geom->gCon[0][i]*geom->alpha;
  }
}

void setBCon(const struct geometry* restrict geom,
             struct fluidElement* restrict elem)
{
  REAL uCov[NDIM];

  conToCov(elem->uCon, geom, uCov);

  elem->bCon[0] =   elem->primVars[B1]*elem->uCov[1]
                  + elem->primVars[B2]*elem->uCov[2] 
                  + elem->primVars[B3]*elem->uCov[3];
  
  for (int i=1; i<NDIM; i++)
  {
    elem->bCon[i] = (  elem->primVars[U3+i] 
                     + elem->bCon[0]*elem->uCon[i]
                    )/elem->uCon[0];
  }
}

void setFluidElement(const REAL primVars[DOF],
                     const struct geometry* restrict geom,
                     struct fluidElement* restrict elem)
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
}

void matterCurrent(const struct fluidElement* restrict elem,
                   const struct geometry* restrict geom,
                   REAL NUp[NDIM])
{
  for (int mu=0; mu<NDIM; mu++)
  {
    NUp[mu] = elem->primVars[RHO]*elem->uCon[mu];
  }
}

void stressTensor(const struct fluidElement* restrict elem,
                  const struct geometry* restrict geom,
                  REAL TUpDown[NDIM][NDIM])
{
  REAL pressure = (ADIABATIC_INDEX - 1.)*elem->primVars[UU];
  REAL bCov[NDIM], uCov[NDIM], bSqr;
  
  conToCov(elem->bCon, geom, bCov);
  bSqr = covDotCon(bCov, elem->bCon);

  conToCov(elem->uCon, geom, uCov);

  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      TUpDown[mu][nu] =   (  elem->primVars[RHO] + elem->primVars[UU] 
                           + pressure + bSqr
                          )*elem->uCon[mu]*uCov[nu]

                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)

                        - elem->bCon[mu]*bCov[nu];
    }
  }

}

void computeFluxesAndConservedVars(const struct fluidElement* restrict elem,
                                   const struct geometry* restrict geom,
                                   const int dir,
                                   REAL fluxes[DOF],
                                   REAL conservedVars[DOF])
{
  REAL NUp[NDIM], TUpDown[NDIM][NDIM];
  REAL g = sqrt(-geom->gDet);

  matterCurrent(elem, geom, NUp);
  stressTensor(elem, geom, TUpDown);
 
  fluxes[RHO] = g*NUp[dir];

  fluxes[UU] = g*TUpDown[dir][0];
  fluxes[U1] = g*TUpDown[dir][1];
  fluxes[U2] = g*TUpDown[dir][2];
  fluxes[U3] = g*TUpDown[dir][3];

  fluxes[B1] = g*(elem->bCon[1]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[1]);
  fluxes[B2] = g*(elem->bCon[2]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[2]);
  fluxes[B3] = g*(elem->bCon[3]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[3]);

  conservedVars[RHO] = g*NUp[0];

  conservedVars[UU] = g*TUpDown[0][0];
  conservedVars[U1] = g*TUpDown[0][1];
  conservedVars[U2] = g*TUpDown[0][2];
  conservedVars[U3] = g*TUpDown[0][3];

  conservedVars[B1] = g*(  elem->bCon[1]*elem->uCon[0] 
                         - elem->bCon[0]*elem->uCon[1]);

  conservedVars[B2] = g*(  elem->bCon[2]*elem->uCon[0] 
                         - elem->bCon[0]*elem->uCon[2]);

  conservedVars[B3] = g*(  elem->bCon[3]*elem->uCon[0] 
                         - elem->bCon[0]*elem->uCon[3]);
}

/*  Returns sqrt(-gDet) * T^kappa_lamda * Gamma^lamda_nu_kappa for each nu (Eqn
 *  (4) of HARM paper).
 * 
 *  = sqrt(-gDet)*T^kappa_lamda * g^lamda^mu * Gamma_mu_nu_kappa
 *
*/

void computeSourceTerms(const struct fluidElement* restrict elem,
                        const struct geometry* restrict geom,
                        const REAL X[NDIM],
                        REAL sourceTerms[DOF])
{
  REAL g = geom->sqrt(-gDet);

  sourceTerms[RHO] = 0.;
  sourceTerms[B1] = 0.;
  sourceTerms[B2] = 0.;
  sourceTerms[B3] = 0.;

  for (int nu=0; nu<NDIM; nu++)
  {
    sourceTerms[RHO+nu] = 0.;
    for (int kappa=0; kappa<NDIM; kappa++)
    {
      for (int lamda=0; lamda<NDIM; lamda++)
      {
        for (int mu=0; mu<NDIM; mu++)
        {
          sourceTerms[RHO+nu] +=   
               g * TUpDown[kappa][lamda] * geom->gCon[lamda][mu]
                 * gammaDownDownDown(mu, nu, kappa, X);
        }
      }
    }
  }

}
