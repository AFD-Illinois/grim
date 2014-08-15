#ifndef GRIM_PROBLEM_H_
#define GRIM_PROBLEM_H_

#include <petsc.h>
#include "../inputs.h"
#include "../timestepper/timestepper.h"
#include "../gridzone/gridzone.h"
#include "../boundary/boundary.h"
#include "../physics/physics.h"

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1]);

void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE]);

void problemDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

#endif /* GRIM_PROBLEM_H_ */
