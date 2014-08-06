#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

#include "../inputs.h"
#include "../geometry/geometry.h"

/* Contains all the variables needed for physics. Independent variables are only
 * primVars. The rest are auxiliary variables stored for convenience.
 */

#define DELTA(mu, nu) (mu==nu ? 1 : 0)

struct fluidElement
{
  REAL gamma;
  REAL primVars[DOF];
  REAL uCon[NDIM];
  REAL bCon[NDIM];
};

/* Public functions: */
void setFluidElement(const REAL primVars[DOF],
                     const struct geometry* restrict geom,
                     struct fluidElement* restrict elem);

void matterCurrent(const struct fluidElement* restrict elem,
                   const struct geometry* restrict geom,
                   REAL NUp[NDIM]);

void stressTensor(const struct fluidElement* restrict elem,
                  const struct geometry* restrict geom,
                  REAL TUpDown[NDIM][NDIM]);

void computeFluxesAndConservedVars(const struct fluidElement* restrict elem,
                                   const struct geometry* restrict geom,
                                   const int dir,
                                   REAL fluxes[DOF],
                                   REAL conservedVars[DOF]);

void computeSourceTerms(const struct fluidElement* restrict elem,
                        const struct geometry* restrict geom,
                        const REAL X[NDIM],
                        REAL sourceTerms[DOF]);

/* Internal functions */
void setGamma(const struct geometry* restrict geom,
              struct fluidElement* restrict elem);

void setUCon(const struct geometry* restrict geom,
             struct fluidElement* restrict elem);

void setBCon(const struct geometry* restrict geom,
             struct fluidElement* restrict elem);

#endif /* GRIM_PHYSICS_H_ */
