#ifndef GRIM_PROBLEM_H_
#define GRIM_PROBLEM_H_

#include <petsc.h>
#include "../inputs.h"
#include "../timestepper/timestepper.h"
#include "../grid/grid.h"
#include "../geometry/geometry.h"
#include "../boundary/boundary.h"
#include "../physics/physics.h"
#include PROBLEM_DATA /* Includes the problem struct defined in the specific
                         problem folder */

//void initialConditions(struct timeStepper ts[ARRAY_ARGS 1]);
//
//void applyFloor(const int iTile, const int jTile,
//                const int X1Start, const int X2Start,
//                const int X1Size, const int X2Size,
//                REAL primTile[ARRAY_ARGS TILE_SIZE]);
//
//void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
//                                     const int X1Start, const int X2Start,
//                                     const int X1Size, const int X2Size,
//                                     REAL primTile[ARRAY_ARGS TILE_SIZE]);
//
//void applyAdditionalProblemSpecificBCs
//(
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  const struct problemData problemSpecificData[ARRAY_ARGS 1],
//  REAL primTile[ARRAY_ARGS TILE_SIZE]
//);
//
//void applyProblemSpecificFluxFilter
//(
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  const struct problemData problemSpecificData[ARRAY_ARGS 1],
//  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
//  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
//);
//
//void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);
//
//void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);
//
//#if (CONDUCTION)
//  void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
//                               struct fluidElement elem[ARRAY_ARGS 1]);
//#endif
//#if (VISCOSITY)
//  void setViscosityParameters(const struct geometry geom[ARRAY_ARGS 1],
//                               struct fluidElement elem[ARRAY_ARGS 1]);
//#endif
//
#endif /* GRIM_PROBLEM_H_ */
