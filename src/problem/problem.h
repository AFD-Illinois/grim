#ifndef GRIM_PROBLEM_H_
#define GRIM_PROBLEM_H_

#include <petsc.h>
#include "../inputs.h"
#include "../timestepper/timestepper.h"
#include "../gridzone/gridzone.h"
#include "../boundary/boundary.h"

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1]);

void applyFloor(REAL primTile[ARRAY_ARGS TILE_SIZE]);

void postStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

#endif /* GRIM_PROBLEM_H_ */
