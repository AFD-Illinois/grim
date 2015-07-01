#include "physics.h"

/* Compute the Lorentz factor \gamma in the coordinate frame
 *
 * gamma = \sqrt{1 + q^2}
 *
 * where
 *
 * q^2 = g_{ij} \~{u}^i \~{u}^j
 *
 * and \~{u}^i are the primitive variables U1, U2, and U3
 *
 * Ref: Eqn (17) and the paragraph below in 
 *      "A Measurement of the Electromagnetic Luminosity of a Kerr Black Hole"
 *      Jonathan C. McKinney and Charles F. Gammie, 2004
 *
 * @param Input: geom, a geometry struct
 * @param Output: elem, a fluid element whose gamma parameter will be set by
 *                this function
*/ 
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

/* Compute the four-velocity u^\mu in the coordinate frame
 *
 * The primitive variables are 
 * U1 \equiv \~{u}^1
 * U2 \equiv \~{u}^2
 * U3 \equiv \~{u}^3
 *
 * with range = (-infty, infty)
 *
 * The variables \~{u}^i are defined using
 *
 * \~{u}^i \equiv u^i + \frac{\gamma \beta^i}{\alpha}
 *
 * where \beta^i = g^{0i} \alpha^2 is the shift
 *       \alpha^2 = -1/g^{00} is the lapse.
 *
 * Therefore, the four-velocity u^i is
 *
 * u^{i} = \~{u}^i - \frac{\gamma \beta^i}{\alpha}
 *       = \~{u}^i - \gamma g^{0i} \alpha
 *
 * In addition
 *
 * u^{0} = \frac{\gamma}{\alpha}
 *
 * Ref: 1) Eqn (17) and the paragraph below in 
 *         "A Measurement of the Electromagnetic Luminosity of a Kerr Black Hole",
 *         Jonathan C. McKinney and Charles F. Gammie, 2004
 *
 *      2) Paragraph above Eqn (16) and above Eqn (22) in 
 *         "Primitive Variable Solvers for Conservative General Relativistic
 *         Magnetohydrodynamics", Noble et. al., 2006
 *
 * @param Input: geom, a geometry struct
 * @param Output: elem, a fluid element whose uCon will be set by this function
 *
*/ 
void setUCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1])
{
  elem->uCon[0] = elem->gamma/geom->alpha;

  for (int i=1; i<NDIM; i++)
  {
    elem->uCon[i] =  elem->primVars[U1+i-1] 
                   - elem->gamma*geom->gCon[0][i]*geom->alpha;
  }
}

/* Compute the magnetic field four-vector b^\mu in the coordinate frame
 *
 * The magnetic field primitive variables are 
 * B1 \equiv *F^{1t}
 * B2 \equiv *F^{2t}
 * B3 \equiv *F^{3t}
 *
 * The variables B^i are used to construct a four-vector b^\mu with
 *
 * b^t = g_{i \mu} B^i u^\mu = B^i u_i
 *
 * b^i = (B^i + u^i b^t)/u^t
 *
 * The above b^\mu is (0, B1, B2, B3) in a comoving tetrad and satisfies 
 * b^\mu u_\mu = 0
 *
 * Ref: 1) Paragraph above Eqn (18) in
 *         "A Measurement of the Electromagnetic Luminosity of a Kerr Black Hole",
 *         Jonathan C. McKinney and Charles F. Gammie, 2004
 *
 * @param Input: geom, a geometry struct
 * @param Output: elem, a fluid element whose bCon will be set by this function
*/ 
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
    elem->bCon[i] = (  elem->primVars[B1+i-1] 
                     + elem->bCon[0]*elem->uCon[i]
                    )/elem->uCon[0];
  }
}

/* Set the {\tt fluidElement} struct, given the primitive variables and the
 * geometry as an input. The fluidElement here represents a fluid element at a
 * particular point in space-time, i.e. it is {\em Eulerian}. Thus every grid
 * zone in the computational domain has a {\tt fluidElement}, as well as a {\tt
 * geometry}. The {\tt fluidElement} struct contains all physics related
 * quantities.
 *
 * @param Input: primVars, an array of all the primitive variables, that we
 *               eventually need to solve for. All physics quantities can be 
 *               recovered by using the primitive variables, along with the 
 *               geometry at a grid zone.
 * @param Input: geom, a geometry struct
 * @param Output: elem, a fluid element whose bCon will be set by this function
*/ 
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
  #if (CONDUCTION || VISCOSITY)
    //setDiffusionCoefficients(geom, elem);
  #endif
  computeMoments(geom, elem);
}

/* Compute and set the various moments of the distribution function in
 * {\tt fluidElement} : 
 * 1) Matter current NUp
 * 2) Stress tensor TUpDown
 *
 * @param Input: geom, a geometry struct
 * @param Output: elem, a fluid element whose  will be set by this function
*/ 
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

/* {\tt grim} solves equations in the following form
 *
 *  dU_i/dt + div.F_i = S_i
 *
 *  i.e.
 *
 *  \partial_\mu F^\mu_i = S_i
 *
 *  where i   = 1..DOF   (iterates over the total number of equations)
 *        \mu = 0,1,2,3  (iterates over all the directions)
 *
 *  This routine computes the spatial fluxes when \mu (given by the input
 *  parameter dir) = 1, 2, 3 and the conserved variables when \mu = 0
 *
 *  Ref: Eqn (2), (4) and (18) in 
 *      "HARM: A Numerical Scheme for General Relativistic Magnetohydrodynamics"
 *      -- Charles F. Gammie, Jonathan C. Mckinney and Gabor Toth, 2003
 *
 *      Also see Eqn (29) and the accompanying paragraph in the above ref for
 *      notes on the modification of the energy equations (fluxes[UU]).
 *
 * @param Input: elem, a {\tt fluidElement} whose moments have been computed.
 * @param Input: geom, a {\tt geometry} struct.
 * @param Input: dir, direction in which the fluxes are needed (sets \mu).
 * @param Output: fluxes, an array containing the fluxes in the dir direction.
*/ 
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

/* {\tt grim} solves equations in the following form
 *
 *  dU_i/dt + div.F_i = S_i
 *
 *  i.e.
 *
 *  \partial_\mu F^\mu_i = S_i
 *
 *  where i   = 1..DOF   (iterates over the total number of equations)
 *        \mu = 0,1,2,3  (iterates over all the directions)
 *
 * This function computes the source terms S_i for all the equations that do NOT
 * have derivatives in them, such as in ideal GRMHD. In this case it returns
 * sqrt(-gDet) * T^kappa^lamda * Gamma_lamda_nu_kappa for each nu. 
 * (Eqn (4) of HARM paper). 
 *
 * When the source terms have spatial and temporal derivatives in them such as
 * in the EMHD model, the source terms are computed using () and added to the
 * residual in computeResidual() in residual.c
 *
 *  Ref:1) Ideal GRMHD source terms: Eqn (4) in 
 *      "HARM: A Numerical Scheme for General Relativistic Magnetohydrodynamics"
 *      -- Charles F. Gammie, Jonathan C. Mckinney and Gabor Toth, 2003
 *
 *      2) Additional source terms from the EMHD model : Eqn () in "grim" 
 *      -- Mani Chandra, Francois Foucart and Charles F. Gammie 
 *
 * @param Input: elem, a {\tt fluidElement} whose moments have been computed.
 * @param Input: geom, a {\tt geometry} struct.
 * @param Input: christoffel, Christoffel symbols of the second kind, indexed
 *               using the GAMMA_UP_DOWN_DOWN macro.
 * @param Output: sources, an array containing the sources terms for all the
 *                equations. Note that only the source terms of equations that
 *                do NOT have spatial and temporal derivatives in the source
 *                terms are computed. The rest are set to zero, and are added
 *                later in computeResidual() in residual.c
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
}

/* Auxiliary function that returns b^\mu b_\mu where b^\mu is defined in
 * setBCon() above.
 *
 * @param Input: elem, a {\tt fluidElement} whose bCon has been set.
 * @param Input: geom, a {\tt geometry} struct.
 * @param Output: bSqr, bSqr = b^\mu b_\mu
 */
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
