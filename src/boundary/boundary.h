#ifndef GRIM_BOUNDARY_H_
#define GRIM_BOUNDARY_H_

#include "../inputs.h"
#include "../gridzone/gridzone.h"
#include "../physics/physics.h"

#define OUTFLOW                       (0)
#define MIRROR                        (1)
#define PERIODIC                      (2)
#define DIRICHLET                     (3)


/* Periodic boundary conditions are handled by Petsc since data needs to be
 * copied if running on more than 1 MPI node. Petsc does all that and puts the
 * appropriate data in the local array. Simply copy that */


/** TILE_BOUNDARY flag set in gridZone2D struct:
 *
 *                      A                             C
 *  +---------+---------+---------+---------+---------+
 *  |         |         |         |         |         | 
 *  |    G2   |    G1   |    P        TILE            |
 *  |         |         |         |         |         | 
 *  +---------+---------+---------+---------+---------+
 *                      B                             D
 *
 *  P       --      current zone
 *  G1      --      Tile ghost zone
 *  G2      --      Tile ghost zone
 *
 *  The entire domain of N1 x N2 physical zones are divided into tiles which are
 *  small enough to fit into the cache of the compute node or of a potential
 *  accelerator. If there are NTile tiles covering the global grid, then the
 *  memory access to RAM is NTile(+ boundaries of the tiles) instead of (N1 x
 *  N2). Each tile has (TILE_SIZE_X1 + 2*NG) x (TILE_SIZE_X2 + 2*NG) number of
 *  grid zones. In particular, the ghost zones during evolution are only for the
 *  tiles. The global grid has no ghost zones. The boundary flags really only
 *  matter for grid zones at the edges of the tiles. A tile can be in the bulk
 *  of the domain or one or more of its edges could coincide with the physical
 *  boundary. If a grid zone at the edge of a tile belongs to a tile whose
 *  boundary does not coincide with a physical boundary, then it is marked as
 *  TILE_BOUNDARY. In this case we need to fill up the ghost zones of the tile
 *  with data from the global grid and this needs a trip to the RAM. On the
 *  other hand if an edge of the tile does indeed coincide with a physical
 *  boundary then the physical boundary flags are set for those particular
 *  edges.
 *
 *  If the edge 'AB' of the current zone 'P' is marked as a TILE_BOUNDARY then
 *  get data for G1 and G2 from the global grid.
 *
 *  OUTFLOW boundary flag set in gridZone2D struct:
 *
 *              A
 *  +-----+-----+-----+-----
 *  |     |     |     |
 *  |  G2 |  G1 |  P  |
 *  |     |     |     |
 *  +-----+-----+-----+-----
 *              B
 *
 *  If the edge 'AB' is marked as OUTFLOW, then copy the variable values from
 *  the physical zone P into the ghost zones G1 and G2.
 *
 *  MIRROR boundary flag set in gridZone2D struct:
 *
 *              A
 *  +-----+-----+-----+-----+
 *  |     |     |     |     |
 *  |  G2 |  G1 |  P1 |  P2 |
 *  |     |     |     |     |
 *  +-----+-----+-----+-----+
 *              B
 * 
 *  If the edge 'AB' is marked as MIRROR, then copy the variable values from the
 *  physical zone P1 into the ghost zone G1, from P2 to G2 and so on.
 *
 *
 *  CONSTANT boundary flag set in gridZone2D struct:
 *
 *              A
 *  +-----+-----+-----+-----
 *  |     |     |     |
 *  |  G2 |  G1 |  P  |
 *  |     |     |     |
 *  +-----+-----+-----+-----
 *              B
 *
 *  If the edge 'AB' is marked as CONSTANT, then set the variable values in the
 *  ghost zones G1 and G2 to a specific values that are provided.
 *
 *  PERIODIC boundary flag set in gridZone2D struct:
 *
 *              A                             C     
 *  +-----+-----+-----+----- _ _ _ -----+-----+-----+-----+
 *  |     |     |     |                 |     |     |     |
 *  |  G2 |  G1 |  Pa |  Pb          Pd |  Pc |  G3 |  G4 |
 *  |     |     |     |      _ _ _      |     |     |     |
 *  +-----+-----+-----+-----       -----+-----+-----+-----+
 *              B                             D
 *
 *  If the edge 'AB' has the same state as 'CD', then copy Pc into G1 and Pd
 *  into G2. Likewise copy Pa into G3 and Pb into G4.
 *
 *  NONE boundary flag set in gridZone2D struct:
 *       
 *         |     |
 *         |     |
 *    -----+-----+-----
 *         |     |
 *         |  P  |
 *         |     |
 *    -----+-----+-----
 *         |     |
 *         |     |
 *
 *  Zone does not share its edge with either the tile boundary or the physical
 *  boundary. It is in the bulk of the domain.
*/

void applyTileBoundaryConditions(const int iTile, const int jTile,
                                 const int X1Start, const int X2Start,
                                 const int X1Size, const int X2Size,
                                 REAL tile[ARRAY_ARGS TILE_SIZE]);
#endif /* GRIM_BOUNDARY_H_ */
