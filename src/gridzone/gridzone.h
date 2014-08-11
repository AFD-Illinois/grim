#ifndef GRIM_GRIDZONE_H_
#define GRIM_GRIDZONE_H_

#include <math.h>
#include "../inputs.h"
#include "../geometry/geometry.h"

#define X1_A (0.)
#define X2_A (0.)
#define X1_B (0.)
#define X2_B (0.)

#define CENTER           (0)
#define FACE_X1          (1)
#define FACE_X2          (2)
#define CORNER           (3)
#define FACE_X1_PLUS_ONE (4)
#define FACE_X2_PLUS_ONE (5)

/* Entire domain is divided into three parts:
 * 1) Global -- The full domain without ghost zones. Reading the global domain
 *              will give the part of the domain that the current processor can
 *              read. Using this means that the associated array does not need
 *              any communication. Used for arrays associated with global Petsc
 *              Vecs.
 *
 * 2) Local -- The part of the full domain that the current processor can read
 *             with ghost zones. Arrays marked local need communication. Used
 *             for arrays associated with local Petsc Vecs.
 *
 * 3) Tile -- Small chunk of the local/global array that can fit inside the
 *            cache. Used for arrays created on our own.
 */
#if (COMPUTE_DIM==1)

#define TILE_SIZE ((TILE_SIZE_X1+2*NG)*DOF)
#define INDEX_TILE(zone,var) (var + DOF*(zone.iInTile + NG) )

#define INDEX_TILE_PLUS_ONE_X1(zone,var) (var + DOF*(zone.iInTile + 1 + NG) )

#define INDEX_TILE_PLUS_TWO_X1(zone,var) (var + DOF*(zone.iInTile + 2 + NG ) )

#define INDEX_TILE_MINUS_ONE_X1(zone,var) (var + DOF*(zone.iInTile - 1 + NG) )

#define INDEX_TILE_MINUS_TWO_X1(zone,var) (var + DOF*(zone.iInTile - 2 + NG) )

#define INDEX_TILE_MANUAL(iInTile,jInTile,var) (var + DOF*(iInTile + NG) )

#define INDEX_LOCAL(zone,var) (var + DOF*(zone.i))
#define INDEX_LOCAL_MANUAL(i,j,zone,var) (var + DOF*(i) )

#define INDEX_GLOBAL(zone,var) (var + DOF*(zone.i))
#define INDEX_GLOBAL_MANUAL(i,j,zone,var) (var + DOF*(i) )


#elif (COMPUTE_DIM==2)

#define TILE_SIZE ((TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*DOF)
#define INDEX_TILE(zone,var) (var + DOF*(\
                                          zone.iInTile + NG \
                                         +(TILE_SIZE_X1+2*NG)\
                                         *(zone.jInTile + NG)\
                                        )\
                             )

#define INDEX_TILE_PLUS_ONE_X1(zone,var) (var + DOF*(\
                                                      zone.iInTile + 1 + NG \
                                                     +(TILE_SIZE_X1+2*NG)\
                                                     *(zone.jInTile + NG)\
                                                    )\
                                         )

#define INDEX_TILE_PLUS_TWO_X1(zone,var) (var + DOF*(\
                                                      zone.iInTile + 2 + NG \
                                                     +(TILE_SIZE_X1+2*NG)\
                                                     *(zone.jInTile + NG)\
                                                    )\
                                         )

#define INDEX_TILE_MINUS_ONE_X1(zone,var) (var + DOF*(\
                                                      zone.iInTile - 1 + NG \
                                                      +(TILE_SIZE_X1+2*NG)\
                                                      *(zone.jInTile + NG)\
                                                     )\
                                          )

#define INDEX_TILE_MINUS_TWO_X1(zone,var) (var + DOF*(\
                                                      zone.iInTile - 2 + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *(zone.jInTile + NG)\
                                                     )\
                                          )

#define INDEX_TILE_PLUS_ONE_X2(zone,var) (var + DOF*(\
                                                      zone.iInTile + NG \
                                                     +(TILE_SIZE_X1 + 2*NG)\
                                                     *(zone.jInTile + 1 + NG)\
                                                    )\
                                         )

#define INDEX_TILE_PLUS_TWO_X2(zone,var) (var + DOF*(\
                                                      zone.iInTile + NG \
                                                     +(TILE_SIZE_X1 + 2*NG)\
                                                     *(zone.jInTile + 2 + NG)\
                                                    )\
                                         )

#define INDEX_TILE_MINUS_ONE_X2(zone,var) (var + DOF*(\
                                                      zone.iInTile + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *(zone.jInTile - 1 + NG)\
                                                     )\
                                          )

#define INDEX_TILE_MINUS_TWO_X2(zone,var) (var + DOF*(\
                                                      zone.iInTile + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *(zone.jInTile - 2 + NG)\
                                                     )\
                                          )

#define INDEX_TILE_MANUAL(iInTile,jInTile,var) \
                                          (var + DOF*(\
                                                      iInTile + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *(jInTile + NG)\
                                                     )\
                                          )

#define INDEX_LOCAL(zone,var) (var + DOF*(zone.i + (\
                                                     zone.X1Size+2*NG\
                                                   )\
                                                  *(zone.j) \
                                         ) \
                              )
#define INDEX_LOCAL_MANUAL(i,j,zone,var)\
                              (var + DOF*(i + (\
                                               zone.X1Size+2*NG\
                                              )\
                                             *(j) \
                                         ) \
                              )

#define INDEX_GLOBAL(zone,var) (var + DOF*(zone.i + zone.X1Size*(zone.j)))
#define INDEX_GLOBAL_MANUAL(i,j,zone,var) (var + DOF*(i + zone.X1Size*(j)))

#endif


struct gridZone
{
  int i, j;
  int iTile, jTile; 
  int iInTile, jInTile;
  int X1Size, X2Size;

  REAL dX1, dX2;

  /* State of the edges: */
  /* OUTFLOW, MIRROR, PERIODIC, CONSTANT */
  int leftEdge, rightEdge;
  int bottomEdge, topEdge;
};

void getXCoords(const struct gridZone zone[ARRAY_ARGS 1],
                const int location,
                REAL X[ARRAY_ARGS NDIM]);

void setGridZone(const int iTile,
                 const int jTile,
                 const int iInTile,
                 const int jIntTile,
                 const int X1Start,
                 const int X2Start,
                 const int X1Size,
                 const int X2Size,
                 struct gridZone zone[ARRAY_ARGS 1]);

#endif /* GRIM_GRIDZONE_H_ */
