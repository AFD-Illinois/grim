#ifndef GRIM_BOUNDARY_H_
#define GRIM_BOUNDARY_H_

#include "../inputs.h"
#include "../grid/grid.h"
#include "../physics/macros.h" /* Determines DOF */
#include "macros.h"

/* Periodic boundary conditions are handled by Petsc since data needs to be
 * copied if running on more than 1 MPI node. Petsc does all that and puts the
 * appropriate data in the local array. Simply copy that */


void applyTileBoundaryConditions(const struct gridTile tile[ARRAY_ARGS 1],
                                 REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)]
                                );
#endif /* GRIM_BOUNDARY_H_ */
