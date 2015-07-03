#ifndef GRIM_GRID_H_
#define GRIM_GRID_H_

#include <math.h>
#include <petsc.h>
#include "../inputs.h"
#include "../boundary/macros.h"
#include "../reconstruct/macros.h"
#include "../physics/macros.h"
#include "macros.h"

struct gridData
{
  int iLocalStart, iLocalSize;
  int jLocalStart, jLocalSize;
  int kLocalStart, kLocalSize;

  int numVar;
  int numTilesX1, numTilesX2;

  DM dm;      Vec vec;
  DM dmGhost; Vec vecGhost;

  POINTER_TO_VEC(ptr);
  POINTER_TO_VEC(ptrGhost);
};

struct gridTile
{
  int iTile, jTile, kGlobal;
  int iLocalStart, jLocalStart, kLocalStart;
};

struct gridZone
{
  int iGlobal, jGlobal, kGlobal;
  int iInTile, jInTile, kInTile;

  REAL dX[NDIM];
};

/* Public functions */
void getXCoords(const struct gridZone zone[ARRAY_ARGS 1],
                const int location,
                REAL X[ARRAY_ARGS NDIM]
               );

void initGridData(const int numVar, const int numGhost,
                  struct gridData grid[ARRAY_ARGS 1]
                 );

void destroyGridData(struct gridData grid[ARRAY_ARGS 1]);

void setPointerToVec(struct gridData grid[ARRAY_ARGS 1]);

void restorePointerToVec(struct gridData grid[ARRAY_ARGS 1]);

void setPointerToVecGhost(struct gridData grid[ARRAY_ARGS 1]);

void restorePointerToVecGhost(struct gridData grid[ARRAY_ARGS 1]);

void setPointerToExternalVec(Vec externalVec,
                             struct gridData grid[ARRAY_ARGS 1]
                            );

void restorePointerToExternalVec(Vec externalVec,
                                 struct gridData grid[ARRAY_ARGS 1]
                                );

void startFillingVecGhost(struct gridData grid[ARRAY_ARGS 1]);

void finishFillingVecGhost(struct gridData grid[ARRAY_ARGS 1]);

void startFillingVecGhostWithExternalVec(Vec externalVec,
                                         struct gridData grid[ARRAY_ARGS 1]
                                        );

void finishFillingVecGhostWithExternalVec(Vec externalVec,
                                          struct gridData grid[ARRAY_ARGS 1]
                                         );

void setTile(const int iTile, const int jTile, const int kGlobal,
             const struct gridData grid[ARRAY_ARGS 1],
             struct gridTile tile[ARRAY_ARGS 1]
            );

void setZone(const int iInTile, const int jInTile, const int kInTile,
             const struct gridTile tile[ARRAY_ARGS 1],
             struct gridZone zone[ARRAY_ARGS 1]
            );

void loadPrimTile(const struct gridData grid[ARRAY_ARGS 1],
                  const struct gridTile tile[ARRAY_ARGS 1],
                  REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)]
                 );

#endif /* GRIM_GRID_H_ */
