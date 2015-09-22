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
  #if (CONDUCTION)
    setConductionParameters(geom, elem);
  #endif
  #if (VISCOSITY)
    setViscosityParameters(geom, elem);
  #endif
  computeMoments(geom, elem);
}

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1])
{
  REAL pressure = (ADIABATIC_INDEX - 1.)*elem->primVars[UU];
  REAL bCov[NDIM], bSqr, uCov[NDIM];
  
  bSqr = getbSqr(elem, geom);

  conToCov(elem->uCon, geom, uCov);
  conToCov(elem->bCon, geom, bCov);

  //Recover dP from psi (which is either dP of dP * (beta/T)^(1/2)
  //depending on whether we use higher order terms in the evolution
  //equation for dP.
  REAL Rho = elem->primVars[RHO];
  if(Rho<RHO_FLOOR_MIN)
    Rho=RHO_FLOOR_MIN;
  REAL U = elem->primVars[UU];
  if(U<UU_FLOOR_MIN)
    U = UU_FLOOR_MIN;
  REAL T = (ADIABATIC_INDEX-1.)*U/Rho;
  if(T<1.e-12)
    T=1.e-12;
  #if  (VISCOSITY)
    REAL dP = elem->primVars[PSI];
    #if (HIGHORDERTERMS_VISCOSITY)
      REAL betaV = elem->tauVis*0.5/elem->eta;
      dP *= sqrt(T/betaV);
    #endif
  #endif

  //Same recovering q from PHI, which is either q or q*(beta/T)^(1/2)
  #if (CONDUCTION)
    REAL q = elem->primVars[PHI];
    #if (HIGHORDERTERMS_CONDUCTION)  
      REAL betaC = elem->tau/elem->kappa/T;
      q *= sqrt(T/betaC);
    #endif
  #endif

  REAL FakeEMHDCoeff = 0.;
  #if (FAKE_EMHD)
    REAL P   = (ADIABATIC_INDEX-1.)*U;
    REAL y = 2.*P/(bSqr+1.e-16);
    FakeEMHDCoeff = 0.5*(-3.*y-2.+sqrt((3.*y+2)*(3.*y+2)+12*y));
  #endif
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
        + q/sqrt(bSqr)
        * (elem->uCon[mu]*bCov[nu] + elem->bCon[mu]*uCov[nu])
      #endif 
      #if  (VISCOSITY)
	//Add -dP*(b^mu b^nu b^{-2} - 1/3 h^{mu nu})
	-dP/bSqr*elem->bCon[mu]*bCov[nu]
	+dP/3.
	*(DELTA(mu, nu)+elem->uCon[mu]*uCov[nu])
      #endif
	-FakeEMHDCoeff*elem->bCon[mu]*bCov[nu]/2.
	+FakeEMHDCoeff*bSqr/6.*(DELTA(mu, nu)+elem->uCon[mu]*uCov[nu])
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

  #if (CONDUCTION)
    fluxes[PHI] = g*(elem->uCon[dir]*elem->primVars[PHI]);
  #endif
  #if (VISCOSITY)
    fluxes[PSI] = g*(elem->uCon[dir]*elem->primVars[PSI]);
  #endif
}

/*  Returns sqrt(-gDet) * T^kappa^lamda * Gamma_lamda_nu_kappa for each nu.
 *  (Eqn (4) of HARM paper)
 * 
 */

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
			const struct gridZone zone[ARRAY_ARGS 1],
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


  #if (ADD_WIND_SOURCE)
    //Add wind source
    REAL XCoords[NDIM];
    getXCoords(zone, CENTER, XCoords);
    REAL xCoords[NDIM];
    XTox(geom->XCoords, xCoords);
    REAL r = xCoords[1];
    REAL cth = cos(xCoords[2]);
    REAL drhodt = 1.e-6*cth*cth*cth*cth/(1.+r*r)/(1+r*r);
    REAL Twind = 10.;
    sourceTerms[RHO]+=drhodt*g/geom->alpha;
    sourceTerms[UU]-=drhodt*g*Twind/(ADIABATIC_INDEX-1.);
  #endif

  #endif
}

REAL getbSqr(const struct fluidElement elem[ARRAY_ARGS 1],
             const struct geometry geom[ARRAY_ARGS 1])
{
  REAL bCov[NDIM], bSqr;
  conToCov(elem->bCon, geom, bCov);
  bSqr = covDotCon(bCov, elem->bCon);
  if (bSqr < 1e-20)
  {
    bSqr = 1e-20;
  }

  return bSqr;
}
