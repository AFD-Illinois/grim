#ifndef GRIM_RECONSTRUCT_H_
#define GRIM_RECONSTRUCT_H_

#include "../inputs.h"
#include "../physics/physics.h"

/* Only Monotonized central slope limiter implemented so far. */
void slopeLim(const REAL left[ARRAY_ARGS DOF],
              const REAL mid[ARRAY_ARGS DOF],
              const REAL right[ARRAY_ARGS DOF],
              REAL ans[ARRAY_ARGS DOF]);

#endif /* GRIM_RECONSTRUCT_H_ */
