//#ifndef GRIM_GRIDZONE_H_
//#define GRIM_GRIDZONE_H_

#include <math.h>
//#include "inputs.h"

/* struct gridZone2D 
 *
 * Data structure holding all the grid information of a specific zone of the
 * global grid, including the location of the zone (i,j) with respect to the
 * global grid, the location of the same zone with respect to the local tile
 * (iTile,jTile) and the size of the zone (dX1, dX2). In addition, the data
 * structure contains boundary information about the zone edges. The four zone
 * edges can be one of either:
 *  1) TILE_BOUNDARY
 *  2) OUTFLOW
 *  3) MIRROR
 *  4) CONSTANT
 *  5) PERIODIC
 *  6) NONE
 * See boundary.h for a detailed description.
*/

#define REAL double
#define NDIM 4
#define CENTER  (10000)
#define FACE_X1 (10001)
#define FACE_X2 (10002)
#define CORNER  (10003)

#define X1_A 0.
#define X2_A 0.
#define H_SLOPE 0.


struct gridZone2D
{
    int i, j;               // Indices of the global domain 0<=i,j<=N1,N2
    int iTile, jTile;       // Indices of the local tile 
                            // 0<=iTile,jTile<=TILE_SIZE_X1, TILE_SIZE_X2

    REAL dX1, dX2;          // dX1 and dX2 of the particular grid zone

    /* State of the edges: */
    /* OUTFLOW, AXISYMMETRIC, PERIODIC, CONSTANT, TILE_BOUNDARY, NONE */
    int leftEdge, rightEdge;
    int bottomEdge, topEdge;
};

void xToX(const REAL x[NDIM], REAL X[NDIM]);

void XTox(const REAL X[NDIM], REAL x[NDIM]);

void getXCoords2D(const struct gridZone2D zone,
                  const int location,
                  REAL X[NDIM]);

//#endif /* GRIM_GRIDZONE_H_ */
