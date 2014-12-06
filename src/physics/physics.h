#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

#include "../inputs.h"
#include "../geometry/geometry.h"

/* Primitive variable mnemonics */
#define RHO             (0)
#define UU              (1)
#define U1              (2)
#define U2              (3)
#define U3              (4)
#define B1              (5)
#define B2              (6)
#define B3              (7)
#if (CONDUCTION)
  #define PHI           (8)
  #define DOF           (9)
#else
  #define DOF           (8)
#endif


/* Contains all the variables needed for physics. Independent variables are only
 * primVars. The rest are auxiliary variables stored for convenience.
 */

#define DELTA(mu, nu) (mu==nu ? 1 : 0)

/* Indices for elem.moments[] */
#define N_UP(mu) (mu)
#define T_UP_DOWN(mu, nu) (nu + NDIM*(mu) + NDIM)

/* Indices for the Christoffel symbols */
#define GAMMA_UP_DOWN_DOWN(eta,mu,nu) (eta+NDIM*(mu+NDIM*(nu) ) )

struct fluidElement
{
  REAL gamma;
  REAL primVars[DOF];
  REAL uCon[NDIM];
  REAL bCon[NDIM];
  REAL moments[20]; 
  /* TODO: We have 4 independent components of N^\mu and 10
  independent components of T^\mu^\nu, so do something later to exploit the
  symmetry of T^\mu^\nu */
};

/* Public functions: */
void setFluidElement(const REAL primVars[ARRAY_ARGS DOF],
                     const struct geometry geom[ARRAY_ARGS 1],
                     struct fluidElement elem[ARRAY_ARGS 1]);

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1]);

void computeFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                   const struct geometry geom[ARRAY_ARGS 1],
                   const int dir,
                   REAL fluxes[ARRAY_ARGS DOF]);

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL christoffel[ARRAY_ARGS 64],
                        REAL sourceTerms[ARRAY_ARGS DOF]);

/* Internal functions */
void setGamma(const struct geometry geom[ARRAY_ARGS 1],
              struct fluidElement elem[ARRAY_ARGS 1]);

void setUCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

void setBCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

#endif /* GRIM_PHYSICS_H_ */
