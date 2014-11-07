#ifndef GRIM_RECONSTRUCT_H_
#define GRIM_RECONSTRUCT_H_

#include "../inputs.h"
#include "../physics/physics.h"

#define X1 (1)
#define X2 (2)

/* Only Monotonized central slope limiter implemented so far. */
void slopeLim(const REAL left[ARRAY_ARGS DOF],
              const REAL mid[ARRAY_ARGS DOF],
              const REAL right[ARRAY_ARGS DOF],
              REAL ans[ARRAY_ARGS DOF]);

void reconstruct(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                 const int dir,
                 const int iTile, const int jTile,
                 const int X1Start, const int X2Start,
                 const int X1Size, const int X2Size,
                 REAL primVarsLeft[ARRAY_ARGS TILE_SIZE],
                 REAL primVarsRight[ARRAY_ARGS TILE_SIZE]);

#endif /* GRIM_RECONSTRUCT_H_ */
