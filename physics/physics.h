#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

/* Contains all the variables needed for physics. Independent variables are only
 * primVars. The rest are auxiliary variables stored for convenience.
 */
struct fluidElement
{
    REAL gamma;
    REAL primVars[DOF];
    REAL uCon[NDIM];
    REAL bCon[NDIM];
};

struct dFluidElementDt
{
    REAL dGammaDt;
    REAL dPrimVarsDt[DOF];
    REAL dUConDt[NDIM];
    REAL dBConDt[NDIM];
}

/* Public functions: */
void setFluidElement(const REAL primVars[DOF],
                     const struct geometry geom,
                     struct fluidElement elem);

void setDFluidElementDt(const struct fluidElement,
                        const REAL dPrimVarsDt[DOF],
                        const struct geometry geom,
                        struct dFluidElementDt dElemDt);

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

void computeDUDt(const struct fluidElement elem,
                 const struct dFluidElementDt dElemDt,
                 const struct geometry geom,
                 REAL dUDt[DOF]);


/* Internal functions */
void setGamma(const struct geometry geom,
              struct fluidElement elem);

void setDGammaDt(const struct geometry geom,
                 const struct fluidElement elem,
                 struct dFluidElementDt dElemDt);

void setUCon(const struct geometry geom,
             struct fluidElement elem);

void setDUConDt(const struct geometry geom,
                const struct fluidElement elem,
                struct dFluidElementDt dElemDt);

void setBCon(const struct geometry geom,
             struct fluidElement elem);

void setDBConDt(const struct geometry geom,
                const struct fluidElement elem,
                struct dFluidElementDt dElemDt);


#endif /* GRIM_PHYSICS_H_ */
