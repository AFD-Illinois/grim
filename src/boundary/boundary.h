#ifndef GRIM_BOUNDARY_H_
#define GRIM_BOUNDARY_H_

#include "../inputs.h"
#include "../grid/grid.h"
#include "../physics/macros.h"
#include "macros.h"

/* Periodic boundary conditions are handled by Petsc since data needs to be
 * copied if running on more than 1 MPI node. Petsc does all that and puts the
 * appropriate data in the local array. Simply copy that */


//void applyTileBoundaryConditions(const int iTile, const int jTile,
//                                 const int X1Start, const int X2Start,
//                                 const int X1Size, const int X2Size,
//                                 REAL tile[ARRAY_ARGS TILE_SIZE]);
#endif /* GRIM_BOUNDARY_H_ */
