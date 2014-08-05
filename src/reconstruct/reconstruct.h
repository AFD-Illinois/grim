#ifndef REAPER_RECONSTRUCT_H_
#define REAPER_RECONSTRUCT_H_

#include "../inputs.h"
#include "../physics/physics.h"

/* Only Monotonized central slope limiter implemented so far. */
void slopeLim(REAL left[DOF],
              REAL mid[DOF],
              REAL right[DOF],
              REAL ans[DOF]);

#endif /* REAPER_RECONSTRUCT_H_ */
