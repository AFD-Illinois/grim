#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

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
                     const struct geometry geom,
                     struct fluidElement elem);

void matterCurrent(const struct fluidElement elem,
                   const struct geometry geom,
                   REAL NUp[NDIM]);

void stressTensor(const struct fluidElement elem,
                  const struct geometry geom,
                  REAL TUpDown[NDIM][NDIM]);

void computeFluxesAndConservedVars(const struct fluidElement elem,
                                   const struct geometry geom,
                                   const int dir,
                                   REAL fluxes[DOF],
                                   REAL conservedVars[DOF]);

/* Internal functions */
void setGamma(const struct geometry geom,
              struct fluidElement elem);

void setUCon(const struct geometry geom,
             struct fluidElement elem);

void setBCon(const struct geometry geom,
             struct fluidElement elem);

#endif /* GRIM_PHYSICS_H_ */
