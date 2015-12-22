#ifndef GRIM_BOUNDARY_H_
#define GRIM_BOUNDARY_H_

#include "../params.hpp"


/* Periodic boundary conditions are handled by Petsc since data needs to be
 * copied if running on more than 1 MPI node. Petsc does all that and puts the
 * appropriate data in the local array. Simply copy that */

#endif /* GRIM_BOUNDARY_H_ */
