#ifndef GRIM_BOUNDARY_H_
#define GRIM_BOUNDARY_H_

#include "../params.hpp"
#include "../grid/grid.hpp"

namespace boundaries
{
  void applyBoundaryConditions(const int boundaryLeft, const int boundaryRight,
                               const int boundaryTop, const int boundaryBottom,
                               const int boundaryFront, const int boundaryBack,
                               grid &prim
                              );
}
/* Periodic boundary conditions are handled by Petsc since data needs to be
 * copied if running on more than 1 MPI node. Petsc does all that and puts the
 * appropriate data in the local array. Simply copy that */

#endif /* GRIM_BOUNDARY_H_ */
