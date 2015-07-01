#ifndef GRIM_RECONSTRUCT_H_
#define GRIM_RECONSTRUCT_H_

#include "../inputs.h"
#include "../physics/macros.h"
#include "macros.h"

/* Contains Monotonized central and MinMod. */
void slopeLim(const REAL left[ARRAY_ARGS DOF],
              const REAL mid[ARRAY_ARGS DOF],
              const REAL right[ARRAY_ARGS DOF],
              REAL ans[ARRAY_ARGS DOF]
             );

void reconstruct(const REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)],
                 const struct gridTile tile[ARRAY_ARGS 1],
                 const int dir,
                 REAL primVarsLeft[ARRAY_ARGS TILE_SIZE(DOF)],
                 REAL primVarsRight[ARRAY_ARGS TILE_SIZE(DOF)]
                );

/* MP5 reconstruction from Suresh and Huynh, 1997. Code taken from PLUTO */
REAL MP5_Reconstruct(const REAL Fjm2,
                     const REAL Fjm1,
                     const REAL Fj,
                     const REAL Fjp1,
                     const REAL Fjp2
                    );

#endif /* GRIM_RECONSTRUCT_H_ */
