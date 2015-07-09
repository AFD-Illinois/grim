#include "physics.h"

void setFluidElement
  (POINTER_TO_VEC(prim),
   const struct gridStrip strip[ARRAY_ARGS 1],
   const struct geometry geom[ARRAY_ARGS 1],
   struct fluidElementStrip elem[ARRAY_ARGS 1]
  )
{
  LOOP_INSIDE_STRIP(0, TILE_SIZE_X1)
  {
    SET_STRIP_GLOBAL_INDICES(strip);

    elem->rho[iInStrip] = INDEX_VEC(prim, strip, RHO);
    elem->uu[iInStrip]  = INDEX_VEC(prim, strip, UU);
    elem->u1[iInStrip]  = INDEX_VEC(prim, strip, U1);
    elem->u2[iInStrip]  = INDEX_VEC(prim, strip, U2);
    elem->u3[iInStrip]  = INDEX_VEC(prim, strip, U3);
    elem->b1[iInStrip]  = INDEX_VEC(prim, strip, B1);
    elem->b2[iInStrip]  = INDEX_VEC(prim, strip, B2);
    elem->b3[iInStrip]  = INDEX_VEC(prim, strip, B3);
  }

  LOOP_INSIDE_STRIP(0, TILE_SIZE_X1)
  {
    elem->gammaLorentz[iInStrip] = 
      sqrt(1 + (geom->gCov11[iInStrip]*elem->u1[iInStrip]*elem->u1[iInStrip])
             + (geom->gCov22[iInStrip]*elem->u2[iInStrip]*elem->u2[iInStrip])
             + (geom->gCov33[iInStrip]*elem->u3[iInStrip]*elem->u3[iInStrip])
             + 2.*(
                    (  geom->gCov12[iInStrip] 
                     * elem->u1[iInStrip] * elem->u2[iInStrip]
                    )
                  + (  geom->gCov13[iInStrip] 
                     * elem->u1[iInStrip] * elem->u3[iInStrip]
                    )
                  + (  geom->gCov23[iInStrip] 
                     * elem->u2[iInStrip] * elem->u3[iInStrip]
                    )
                  )
          );
               
    elem->uUp0[iInStrip] =  elem->gammaLorentz[iInStrip]
                          / geom->alpha[iInStrip];

    elem->uUp1[iInStrip] =  elem->u1[iInStrip]
                          - (  elem->gammaLorentz[iInStrip]
                             * geom->gCon01[iInStrip]
                             * geom->alpha[iInStrip]
                            );

    elem->uUp2[iInStrip] =  elem->u2[iInStrip]
                          - (  elem->gammaLorentz[iInStrip]
                             * geom->gCon02[iInStrip]
                             * geom->alpha[iInStrip]
                            );

    elem->uUp3[iInStrip] =  elem->u3[iInStrip]
                          - (  elem->gammaLorentz[iInStrip]
                             * geom->gCon03[iInStrip]
                             * geom->alpha[iInStrip]
                            );
  }

  LOOP_INSIDE_STRIP(0, TILE_SIZE_X1)
  {
    elem->uDown0[iInStrip] =
            geom->gCov00[iInStrip] * elem->uUp0[iInStrip]
          + geom->gCov01[iInStrip] * elem->uUp1[iInStrip]
          + geom->gCov02[iInStrip] * elem->uUp2[iInStrip]
          + geom->gCov03[iInStrip] * elem->uUp3[iInStrip];

    elem->uDown1[iInStrip] =
            geom->gCov10[iInStrip] * elem->uUp0[iInStrip]
          + geom->gCov11[iInStrip] * elem->uUp1[iInStrip]
          + geom->gCov12[iInStrip] * elem->uUp2[iInStrip]
          + geom->gCov13[iInStrip] * elem->uUp3[iInStrip];

    elem->uDown2[iInStrip] =
            geom->gCov20[iInStrip] * elem->uUp0[iInStrip]
          + geom->gCov21[iInStrip] * elem->uUp1[iInStrip]
          + geom->gCov22[iInStrip] * elem->uUp2[iInStrip]
          + geom->gCov23[iInStrip] * elem->uUp3[iInStrip];

    elem->uDown3[iInStrip] =
            geom->gCov30[iInStrip] * elem->uUp0[iInStrip]
          + geom->gCov31[iInStrip] * elem->uUp1[iInStrip]
          + geom->gCov32[iInStrip] * elem->uUp2[iInStrip]
          + geom->gCov33[iInStrip] * elem->uUp3[iInStrip];
  }

  LOOP_INSIDE_STRIP(0, TILE_SIZE_X1)
  {
    elem->bUp0[iInStrip] =
            elem->b1[iInStrip] * elem->uDown1[iInStrip]
          + elem->b2[iInStrip] * elem->uDown2[iInStrip]
          + elem->b3[iInStrip] * elem->uDown3[iInStrip];

    elem->bUp1[iInStrip] =
        (elem->b1[iInStrip] + elem->bUp0[iInStrip] * elem->uUp1[iInStrip])
      / elem->uUp0[iInStrip];

    elem->bUp2[iInStrip] =
        (elem->b2[iInStrip] + elem->bUp0[iInStrip] * elem->uUp2[iInStrip])
      / elem->uUp0[iInStrip];

    elem->bUp3[iInStrip] =
        (elem->b3[iInStrip] + elem->bUp0[iInStrip] * elem->uUp3[iInStrip])
      / elem->uUp0[iInStrip];
  }

  LOOP_INSIDE_STRIP(0, TILE_SIZE_X1)
  {
    elem->bDown0[iInStrip] =
            geom->gCov00[iInStrip] * elem->bUp0[iInStrip]
          + geom->gCov01[iInStrip] * elem->bUp1[iInStrip]
          + geom->gCov02[iInStrip] * elem->bUp2[iInStrip]
          + geom->gCov03[iInStrip] * elem->bUp3[iInStrip];

    elem->bDown1[iInStrip] =
            geom->gCov10[iInStrip] * elem->bUp0[iInStrip]
          + geom->gCov11[iInStrip] * elem->bUp1[iInStrip]
          + geom->gCov12[iInStrip] * elem->bUp2[iInStrip]
          + geom->gCov13[iInStrip] * elem->bUp3[iInStrip];

    elem->bDown2[iInStrip] =
            geom->gCov20[iInStrip] * elem->bUp0[iInStrip]
          + geom->gCov21[iInStrip] * elem->bUp1[iInStrip]
          + geom->gCov22[iInStrip] * elem->bUp2[iInStrip]
          + geom->gCov23[iInStrip] * elem->bUp3[iInStrip];

    elem->bDown3[iInStrip] =
            geom->gCov30[iInStrip] * elem->bUp0[iInStrip]
          + geom->gCov31[iInStrip] * elem->bUp1[iInStrip]
          + geom->gCov32[iInStrip] * elem->bUp2[iInStrip]
          + geom->gCov33[iInStrip] * elem->bUp3[iInStrip];
  }

  LOOP_INSIDE_STRIP(0, TILE_SIZE_X1)
  {
    elem->bSqr[iInStrip] = 
          elem->bUp0[iInStrip] * elem->bDown0[iInStrip]
        + elem->bUp1[iInStrip] * elem->bDown1[iInStrip]
        + elem->bUp2[iInStrip] * elem->bDown2[iInStrip]
        + elem->bUp3[iInStrip] * elem->bDown3[iInStrip];
  }

  setParameters(elem);
}

void computeConserved(const struct fluidElementStrip elem[ARRAY_ARGS 1],
                      const struct geometryStrip geom[ARRAY_ARGS 1],
                      POINTER_TO_VEC(conserved)
                     )
{
  STATIC_ARRAY(conservedRHO, TILE_SIZE_X1);
  STATIC_ARRAY(conservedUU,  TILE_SIZE_X1);
  STATIC_ARRAY(conservedU1,  TILE_SIZE_X1);
  STATIC_ARRAY(conservedU2,  TILE_SIZE_X1);
  STATIC_ARRAY(conservedU3,  TILE_SIZE_X1);
  STATIC_ARRAY(conservedB1,  TILE_SIZE_X1);
  STATIC_ARRAY(conservedB2,  TILE_SIZE_X1);
  STATIC_ARRAY(conservedB3,  TILE_SIZE_X1);

  LOOP_INSIDE_STRIP(0, TILE_SIZE_X1)
  {
    conservedRHO[INDEX_STRIP(iInStrip)] =  
                      geom->g[INDEX_STRIP(iInStrip)]
                    * elem->rho[INDEX_STRIP(iInStrip)]
                    * elem->uUp0[INDEX_STRIP(iInStrip)];

    conservedUU[INDEX_STRIP(iInStrip)] = 
                      

  }

}

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1]
                   )
{
  REAL pressure = (ADIABATIC_INDEX - 1.)*elem->primVars[UU];
  REAL bCov[NDIM], bSqr, uCov[NDIM];
  
  bSqr = getbSqr(elem, geom);

  conToCov(elem->uCon, geom, uCov);
  conToCov(elem->bCon, geom, bCov);

  /* Heat flux is q. However the variable evolved to solve for the heat flux is
   * PHI.
   * Now the definition of PHI is:
   *      1) q              (when HIGH_ORDER_TERMS_CONDUCTION == OFF)
   * (or) 2) q * sqrt(b1/T) (when HIGH_ORDER_TERMS_CONDUCTION == ON) 
   *
   * where b1 = \tau/(2 \rho \chi) */

  #if  (CONDUCTION)

    #if (HIGH_ORDER_TERMS_CONDUCTION == ON)
      REAL b1 = elem->tauConduction / (2. * elem->primVars[RHO] * elem->chi);
      REAL T = (ADIABATIC_INDEX-1.)*elem->primVars[UU]/elem->primVars[RHO];
      if (T < 1.e-12) T=1.e-12;
      REAL q = elem->primVars[PHI] * sqrt(T/b1);
    #else
      REAL q = elem->primVars[PHI];
    #endif

  #endif


  /* Pressure anisotropy is deltaP. However the variable evolved to solve for
   * pressure anisotropy is PSI.
   * Now the definition of PSI is:
   *      1) deltaP              (when HIGH_ORDER_TERMS_VISOCITY == OFF)
   * (or) 2) deltaP * sqrt(b2/T) (when HIGH_ORDER_TERMS_VISOCITY == ON) 
   *
   * where b2 = \tau/(2 \rho \nu) */

  #if  (VISCOSITY)

    #if (HIGH_ORDER_TERMS_VISCOSITY == ON)
      REAL b2 = elem->tauViscosity / (2. * elem->primVars[RHO] * elem->nu);
      REAL T = (ADIABATIC_INDEX-1.) * elem->primVars[UU] / elem->primVars[RHO];
      if (T < 1.e-12) T = 1.e-12;
      REAL deltaP = elem->primVars[PSI] * sqrt(T/b2);
    #else
      REAL deltaP = elem->primVars[PSI];
    #endif

  #endif

  #pragma vector aligned
  for (int mu=0; mu<NDIM; mu++)
  {
    elem->moments[N_UP(mu)] = elem->primVars[RHO]*elem->uCon[mu];

    #pragma vector aligned
    for (int nu=0; nu<NDIM; nu++)
    {
      elem->moments[T_UP_DOWN(mu,nu)] =   
                          (  elem->primVars[RHO] + elem->primVars[UU]
                           + pressure + bSqr
                          )*elem->uCon[mu]*uCov[nu]

                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)

                        - elem->bCon[mu]*bCov[nu]
      #if (CONDUCTION) 
        /* Add the conduction part of the stress tensor
         *   u^\mu q_\nu + u_\nu q^\mu
         * = u^\mu q \hat{b}_\nu + u_\nu q \hat{b}^\mu */
        + q/sqrt(bSqr)
        * (elem->uCon[mu]*bCov[nu] + elem->bCon[mu]*uCov[nu])
      #endif 
      #if  (VISCOSITY)
	      /* Add the shear stress 
         * \tau^\mu_\nu = -deltaP*(\hat{b}^\mu \hat{b}_\nu - 1/3 h^\mu_\nu)
         *
         * where h^\mu_\nu = \delta^\mu_\nu + u^\mu u_\nu */
	      - deltaP/bSqr * elem->bCon[mu]*bCov[nu]
	      + deltaP/3.   * (DELTA(mu, nu) + elem->uCon[mu]*uCov[nu])
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

  #if (CONDUCTION)
    fluxes[PHI] = g*(elem->uCon[dir]*elem->primVars[PHI]);
  #endif
  #if (VISCOSITY)
    fluxes[PSI] = g*(elem->uCon[dir]*elem->primVars[PSI]);
  #endif
}


void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL christoffel[ARRAY_ARGS 64],
                        REAL sourceTerms[ARRAY_ARGS DOF])
{
  memset(sourceTerms, 0., sizeof(REAL[DOF]));

  #if (METRIC==MODIFIED_KERR_SCHILD)
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
}

REAL getbSqr(const struct fluidElement elem[ARRAY_ARGS 1],
             const struct geometry geom[ARRAY_ARGS 1])
{
  REAL bCov[NDIM], bSqr;
  conToCov(elem->bCon, geom, bCov);
  bSqr = covDotCon(bCov, elem->bCon);
  if (bSqr < 1e-25)
  {
    bSqr = 1e-25;
  }

  return bSqr;
}
