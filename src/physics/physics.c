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
  elem->kappa = 0.1;
  elem->tau = 10.;
}

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1])
{
  REAL pressure = (ADIABATIC_INDEX - 1.)*elem->primVars[UU];
  REAL bCov[NDIM], bSqr, uCov[NDIM];
  
  conToCov(elem->bCon, geom, bCov);
  bSqr = covDotCon(bCov, elem->bCon);
  if (bSqr < 1e-15)
  {
    bSqr = 1e-15;
  }

  conToCov(elem->uCon, geom, uCov);

  for (int mu=0; mu<NDIM; mu++)
  {
    elem->moments[N_UP(mu)] = elem->primVars[RHO]*elem->uCon[mu];

    for (int nu=0; nu<NDIM; nu++)
    {
      elem->moments[T_UP_DOWN(mu,nu)] =   
                          (  elem->primVars[RHO] + elem->primVars[UU]
                           + pressure + bSqr
                          )*elem->uCon[mu]*uCov[nu]

                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)

                        - elem->bCon[mu]*bCov[nu]
      #if (CONDUCTION) 
        + elem->uCon[mu] * elem->primVars[PHI] * bCov[nu]/sqrt(bSqr)
        + uCov[nu] * elem->primVars[PHI] * elem->bCon[mu]/sqrt(bSqr)
      #endif                  
                        ;

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

  fluxes[UU] = g*elem->moments[T_UP_DOWN(dir, 0)] + fluxes[RHO];
  fluxes[U1] = g*elem->moments[T_UP_DOWN(dir, 1)];
  fluxes[U2] = g*elem->moments[T_UP_DOWN(dir, 2)];
  fluxes[U3] = g*elem->moments[T_UP_DOWN(dir, 3)];

  fluxes[B1] = g*(elem->bCon[1]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[1]);
  fluxes[B2] = g*(elem->bCon[2]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[2]);
  fluxes[B3] = g*(elem->bCon[3]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[3]);
}

#if (CONDUCTION)
void computeConductionFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                             const struct geometry geom[ARRAY_ARGS 1],
                             const int dir,
                             REAL fluxes[ARRAY_ARGS DOF])
{
  REAL temperature = (ADIABATIC_INDEX-1.)*elem->primVars[UU]/elem->primVars[RHO];

  for (int var=0; var<DOF; var++)
  {
    fluxes[var] = 0.;
  }

  fluxes[0] = elem->primVars[PHI];
  fluxes[1] = temperature;
  fluxes[2] = elem->uCon[0];
  fluxes[3] = elem->uCon[1];
  fluxes[4] = elem->uCon[2];
  fluxes[5] = elem->uCon[3];
}

void computeConductionFluxesCoefficients(
                             const struct fluidElement elem[ARRAY_ARGS 1],
                             const struct geometry geom[ARRAY_ARGS 1],
                             const int dir,
                             REAL coefficients[ARRAY_ARGS DOF])
{
  REAL bSqr, bCov[NDIM];
  conToCov(elem->bCon, geom, bCov);
  bSqr = covDotCon(bCov, elem->bCon);
  if (bSqr < 1e-15)
  {
    bSqr = 1e-15;
  }
  REAL temperature =   (ADIABATIC_INDEX-1.)
                     * elem->primVars[UU]/elem->primVars[RHO];

  for (int var=0; var<DOF; var++)
  {
    coefficients[var] = 0.;
  }

  coefficients[0] = elem->uCon[dir];
  coefficients[1] = elem->kappa/elem->tau * elem->bCon[dir]/sqrt(bSqr);

  coefficients[2] =   elem->kappa/elem->tau * temperature 
                    * bCov[0]/sqrt(bSqr) * elem->uCon[dir];
  
  coefficients[3] =   elem->kappa/elem->tau * temperature 
                    * bCov[1]/sqrt(bSqr) * elem->uCon[dir];
  
  coefficients[4] =   elem->kappa/elem->tau * temperature 
                    * bCov[2]/sqrt(bSqr) * elem->uCon[dir];

  coefficients[5] =   elem->kappa/elem->tau * temperature 
                    * bCov[3]/sqrt(bSqr) * elem->uCon[dir];

}
#endif

/*  Returns sqrt(-gDet) * T^kappa^lamda * Gamma_lamda_nu_kappa for each nu.
 *  (Eqn (4) of HARM paper)
 * 
*/

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL christoffel[ARRAY_ARGS 64],
                        REAL sourceTerms[ARRAY_ARGS DOF])
{
  for (int var=0; var<DOF; var++)
  {
    sourceTerms[var] = 0.;
  }

  #if (METRIC==KERRSCHILD)
    REAL g = sqrt(-geom->gDet);

    for (int nu=0; nu<NDIM; nu++)
    {
      for (int kappa=0; kappa<NDIM; kappa++)
      {
        for (int lamda=0; lamda<NDIM; lamda++)
        {
          sourceTerms[UU+nu] +=   
                g * elem->moments[T_UP_DOWN(kappa, lamda)]
                  * christoffel[GAMMA_UP_DOWN_DOWN(lamda, kappa, nu)];
        }
      }
    }
  #endif

  #if (CONDUCTION)
    sourceTerms[PHI] = -elem->primVars[PHI]/elem->tau;
    REAL bSqr, bCov[NDIM];
    conToCov(elem->bCon, geom, bCov);
    bSqr = covDotCon(bCov, elem->bCon);
    if (bSqr < 1e-15)
    {
      bSqr = 1e-15;
    }
    REAL temperature =   (ADIABATIC_INDEX-1.)
                       * elem->primVars[UU]/elem->primVars[RHO];

    #if (METRIC==KERRSCHILD)
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int kappa=0; kappa<NDIM; kappa++)
      {
        for (int lamda=0; lamda<NDIM; lamda++)
        {
          sourceTerms[PHI] +=   
           -elem->kappa/elem->tau * temperature * bCov[lamda]/sqrt(bSqr)
          * christoffel[GAMMA_UP_DOWN_DOWN(lamda, kappa, nu)]
          * elem->uCon[kappa] * elem->uCon[nu];
        }
      }
    }

    #endif

  #endif
}
