#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

#include "../inputs.h"
#include "../geometry/geometry.h"   /* Determines METRIC and imports geometry struct*/
#include "../reconstruct/macros.h"  /* Determines NG */
#include "../problem/problem.h"     /* imports setDiffusionCoefficients() */
#include "macros.h"

struct fluidElement
{
  STATIC_ARRAY(rho, TILE_SIZE_X1);
  STATIC_ARRAY(uu,  TILE_SIZE_X1);
  STATIC_ARRAY(u1,  TILE_SIZE_X1);
  STATIC_ARRAY(u2,  TILE_SIZE_X1);
  STATIC_ARRAY(u3,  TILE_SIZE_X1);
  STATIC_ARRAY(b1,  TILE_SIZE_X1);
  STATIC_ARRAY(b2,  TILE_SIZE_X1);
  STATIC_ARRAY(b3,  TILE_SIZE_X1);

  STATIC_ARRAY(gammaLorentz,  TILE_SIZE_X1);
  STATIC_ARRAY(uUp0, TILE_SIZE_X1);
  STATIC_ARRAY(uUp1, TILE_SIZE_X1);
  STATIC_ARRAY(uUp2, TILE_SIZE_X1);
  STATIC_ARRAY(uUp3, TILE_SIZE_X1);

  STATIC_ARRAY(uDown0, TILE_SIZE_X1);
  STATIC_ARRAY(uDown1, TILE_SIZE_X1);
  STATIC_ARRAY(uDown2, TILE_SIZE_X1);
  STATIC_ARRAY(uDown3, TILE_SIZE_X1);

  STATIC_ARRAY(bUp0, TILE_SIZE_X1);
  STATIC_ARRAY(bUp1, TILE_SIZE_X1);
  STATIC_ARRAY(bUp2, TILE_SIZE_X1);
  STATIC_ARRAY(bUp3, TILE_SIZE_X1);

  STATIC_ARRAY(bDown0, TILE_SIZE_X1);
  STATIC_ARRAY(bDown1, TILE_SIZE_X1);
  STATIC_ARRAY(bDown2, TILE_SIZE_X1);
  STATIC_ARRAY(bDown3, TILE_SIZE_X1);
  STATIC_ARRAY(bSqr,   TILE_SIZE_X1);

  STATIC_ARRAY(pressure, TILE_SIZE_X1);

  /* Dissipation coefficients: \chi for conduction and \nu for viscosity
   * Relaxation time scales  : tauConduction and tauViscosity 
   * Ref: EMHD model paper -- Chandra et. al., 2015 */
  #if (CONDUCTION)
    STATIC_ARRAY(chi,           TILE_SIZE_X1);
    STATIC_ARRAY(tauConduction, TILE_SIZE_X1);
  #endif
  #if (VISCOSITY)
    STATIC_ARRAY(nu,           TILE_SIZE_X1);
    STATIC_ARRAY(tauViscosity, TILE_SIZE_X1);
  #endif
};

/* Public functions: */
void setFluidElement
  (const REAL primVars[ARRAY_ARGS DOF][TILE_SIZE_X1],
   const struct geometry geom[ARRAY_ARGS 1],
   struct fluidElement elem[ARRAY_ARGS 1]
  );

void computeFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                   const struct geometry geom[ARRAY_ARGS 1],
                   const int dir,
                   REAL fluxes[ARRAY_ARGS DOF]);

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL christoffel[ARRAY_ARGS 64],
                        REAL sourceTerms[ARRAY_ARGS DOF]);

REAL riemannSolver(const REAL fluxLeft[ARRAY_ARGS DOF],
                   const REAL fluxRight[ARRAY_ARGS DOF],
                   const REAL conservedVarsLeft[ARRAY_ARGS DOF],
                   const REAL conservedVarsRight[ARRAY_ARGS DOF],
                   const REAL primVarsLeft[ARRAY_ARGS DOF],
                   const REAL primVarsRight[ARRAY_ARGS DOF],
                   const struct geometry geom[ARRAY_ARGS DOF],
                   const int dir, REAL fluxes[ARRAY_ARGS DOF]
                  );

void waveSpeeds(const struct fluidElement elem[ARRAY_ARGS 1],
                const struct geometry geom[ARRAY_ARGS 1],
                const int dir,
                REAL cMin[ARRAY_ARGS 1], REAL cMax[ARRAY_ARGS 1]
               );

//#if (CONDUCTION)
//void addConductionSourceTermsToResidual
//(
//  const REAL primTile[ARRAY_ARGS TILE_SIZE],
//  ARRAY(primGlobal), ARRAY(primHalfStepGlobal), ARRAY(primOldGlobal),
//  ARRAY(connectionGlobal),
//  ARRAY(gradTGlobal), ARRAY(graduConGlobal), 
//  ARRAY(graduConHigherOrderTerm1Global),
//  ARRAY(graduConHigherOrderTerm2Global),
//  REAL dt,
//  int computeOldSourceTermsAndOldDivOfFluxes,
//  int computeDivOfFluxAtTimeN,
//  int computeDivOfFluxAtTimeNPlusHalf,
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  ARRAY(residualGlobal)
//);
//
//void computeConductionSpatialGradientTerms
//(
//  const REAL primTile[ARRAY_ARGS TILE_SIZE],
//  const int iTile, const int jTile,
//  const int iInTile, const int jInTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  REAL gradT[COMPUTE_DIM],
//  REAL graduCon[COMPUTE_DIM*NDIM],
//  REAL graduConHigherOrderTerm1[COMPUTE_DIM],
//  REAL graduConHigherOrderTerm2[COMPUTE_DIM]
//);
//#endif
//#if (VISCOSITY)
//void addViscositySourceTermsToResidual
//(
//  const REAL primTile[ARRAY_ARGS TILE_SIZE],
//  ARRAY(primGlobal), ARRAY(primHalfStepGlobal), ARRAY(primOldGlobal),
//  ARRAY(connectionGlobal),
//  ARRAY(graduConVisGlobal), 
//  ARRAY(graduConHigherOrderTerm1VisGlobal),
//  REAL dt,
//  int computeOldSourceTermsAndOldDivOfFluxes,
//  int computeDivOfFluxAtTimeN,
//  int computeDivOfFluxAtTimeNPlusHalf,
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  ARRAY(residualGlobal)
//);
//
//void computeViscositySpatialGradientTerms
//(
//  const REAL primTile[ARRAY_ARGS TILE_SIZE],
//  const int iTile, const int jTile,
//  const int iInTile, const int jInTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  REAL graduConVis[COMPUTE_DIM*NDIM],
//  REAL graduConHigherOrderTerm1Vis[COMPUTE_DIM]
//);
//#endif
//
/* Internal functions */
inline void setGamma(const struct geometry geom[ARRAY_ARGS 1],
              struct fluidElement elem[ARRAY_ARGS 1]);

void setUCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

void setBCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

REAL getbSqr(const struct fluidElement elem[ARRAY_ARGS 1],
             const struct geometry geom[ARRAY_ARGS 1]);
#endif /* GRIM_PHYSICS_H_ */
