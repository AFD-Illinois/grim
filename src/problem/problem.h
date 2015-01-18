#ifndef GRIM_PROBLEM_H_
#define GRIM_PROBLEM_H_

#include <petsc.h>
#include "../inputs.h"
#include "../timestepper/timestepper.h"
#include "../gridzone/gridzone.h"
#include "../geometry/geometry.h"
#include "../boundary/boundary.h"
#include "../physics/physics.h"

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1]);

void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE]);

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
                                     const int X1Start, const int X2Start,
                                     const int X1Size, const int X2Size,
                                     REAL primTile[ARRAY_ARGS TILE_SIZE]);

void applyAdditionalProblemSpecificBCs(const int iTile, const int jTile,
                                       const int X1Start, const int X2Start,
                                       const int X1Size, const int X2Size,
                                       REAL tile[ARRAY_ARGS TILE_SIZE]);

void applyProblemSpecificFluxFilter(const int iTile, const int jTile,
                                    const int X1Start, const int X2Start,
                                    const int X1Size, const int X2Size,
                                    REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
                                    REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]);

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);
void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

#if (CONDUCTION)
void setConductionParameters(struct fluidElement elem[ARRAY_ARGS 1]);
#endif

#endif /* GRIM_PROBLEM_H_ */
