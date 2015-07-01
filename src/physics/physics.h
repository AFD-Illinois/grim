#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

#include "../inputs.h"
#include "../geometry/macros.h"     /* Determines METRIC */
#include "../reconstruct/macros.h"  /* Determines NG */
#include "../problem/problem.h"     /* imports setDiffusionCoefficients() */
#include "macros.h"

/* struct containing all the variables needed for physics. Independent variables
 * are only primVars. The rest are auxiliary variables stored for convenience.
 */
struct fluidElement
{
  REAL primVars[DOF]; /* All the variables that need to be solved for */
  REAL gamma;       /* Lorentz factor      : see setGamma() for description */
  REAL uCon[NDIM];  /* Fluid four-velocity : see setUCon() */
  REAL bCon[NDIM];  /* Magnetic field four-vector : see setBCon() */
  REAL moments[20]; /* 4 components of N^\mu + 16 components of T^{\mu \nu} 
                       : see computeMoments() */

  /* TODO: We have 4 independent components of N^\mu and 10
  independent components of T^\mu^\nu, so do something later to exploit the
  symmetry of T^\mu^\nu.*/

  /* Dissipation coefficients: \chi for conduction and \nu for viscosity
   * Relaxation time scales  : tauConduction and tauViscosity 
   * Ref: EMHD model paper -- Chandra et. al., 2015 */
  #if (CONDUCTION)
    REAL chi, tauConduction;
  #endif
  #if (VISCOSITY)
    REAL nu, tauViscosity;
  #endif
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
void setGamma(const struct geometry geom[ARRAY_ARGS 1],
              struct fluidElement elem[ARRAY_ARGS 1]);

void setUCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

void setBCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

REAL getbSqr(const struct fluidElement elem[ARRAY_ARGS 1],
             const struct geometry geom[ARRAY_ARGS 1]);
#endif /* GRIM_PHYSICS_H_ */
