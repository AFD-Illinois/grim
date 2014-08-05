#include "gridzone.h"

/* Get the coordinates corresponding to a location inside a grid zone.
 *
 *         +---> X1
 *         |
 *         |
 *         v
 *         X2
 *         
 *         (i,j)  (i+0.5,j)
 *           +-----F1-----+
 *           |            |
 *           |            |
 * (i,j+0.5) F2    C      |
 *           |            |
 *           |            |
 *           +------------+
 *
 *      C   --  Center (CENTER)
 *      F1  --  Face in the X1 direction (FACE_X1)
 *      F2  --  Face in the X2 direction (FACE_X2)
 *
 * @param input: zone, Zone in which the coordinates are needed
 * @param input: location, Location inside the zone where the coordinates are needed
 * @param output: X, The coordinates X^mu
*/
void getXCoords(const struct gridZone* restrict zone,
                const int location,
                REAL X[NDIM])
{
#if (COMPUTE_DIM==2)
  if (location==CENTER) // (i+0.5,j+0.5)
  {
    X[0] = 0.; // Don't care about the time coordinate. Just set to 0.
    X[1] = X1_A + (zone->i + 0.5)*zone->dX1;
    X[2] = X2_A + (zone->j + 0.5)*zone->dX2;
    X[3] = 0.;
  } 
  else if (location==FACE_X1) // (i+0.5,j)
  {
    X[0] = 0.;
    X[1] = X1_A + (zone->i + 0.5)*zone->dX1;
    X[2] = X2_A + zone->j*zone->dX2;
    X[3] = 0.;
  }
  else if (location==FACE_X2) // (i,j+0.5)
  {
    X[0] = 0.;
    X[1] = X1_A + zone->i*zone->dX1;
    X[2] = X2_A + (zone->j + 0.5)*zone->dX2;
    X[3] = 0.;
  }
  else if (location==CORNER) // (i,j)
  {
    X[0] = 0.;
    X[1] = X1_A + zone->i*zone->dX1;
    X[2] = X2_A + zone->j*zone->dX2;
    X[3] = 0.;
  }
#elif (COMPUTE_DIM==1)
  if (location==CENTER)
  {
    X[0] = 0.;
    X[1] = X1_A + (zone.i + 0.5)*zone.dX1;
    X[2] = 0.;
    X[3] = 0.;
  }
  else if (location==FACE_X1)
  {
    X[0] = 0.;
    X[1] = X1_A + zone.i*zone.dX1;
    X[2] = 0.;
    X[3] = 0.;
  }
#endif
}

void setGridZone(const int iTile,
                 const int jTile,
                 const int iInTile,
                 const int jIntTile,
                 const int x1Start,
                 const int x2Start,
                 const int x1Size,
                 const int x2Size,
                 struct gridZone* restrict zone)
{
  zone->iTile = iTile;
  zone->jTile = jTile;

  zone->iInTile = iInTile;
  zone->jInTile = jInTile;

  zone->x1Size = x1Size;
  zone->x2Size = x2Size;

  zone->i = x1Start + iInTile + iTile*TILE_SIZE_X1;
  zone->j = x2Start + jIntTile + jTIle*TILE_SIZE_X2;

  zone->dX1 = (X1_B - X1_A)/((REAL)N1);
  zone->dX2 = (X2_B - X2_A)/((REAL)N2);
}
