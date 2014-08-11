#ifndef REAPER_RECONSTRUCT_H_
#define REAPER_RECONSTRUCT_H_

#include "../inputs.h"
#include "../physics/physics.h"

/* Only Monotonized central slope limiter implemented so far. */
void slopeLim(const REAL left[ARRAY_ARGS DOF],
              const REAL mid[ARRAY_ARGS DOF],
              const REAL right[ARRAY_ARGS DOF],
              REAL ans[DOF]);

#endif /* REAPER_RECONSTRUCT_H_ */
