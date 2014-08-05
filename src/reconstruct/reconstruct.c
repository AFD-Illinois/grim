#include "reconstruct.h"

void slopeLim(REAL left[DOF],
              REAL mid[DOF],
              REAL right[DOF],
              REAL ans[DOF])
{
#pragma ivdep
  for (int var=0; var<DOF; var++)
  {
    /* Monotonized Central slope limiter */
    REAL Dqm = 2. * (mid[var] - left[var]);
	  REAL Dqp = 2. * (right[var] - mid[var]);
	  REAL Dqc = 0.5 * (right[var] - left[var]);
	  REAL s = Dqm * Dqp;
	  if (s <= 0.) 
    {
		  ans[var] = 0.;
    }
	  else
    {
		  if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      {
			  ans[var] = Dqm;
      }
		  else if (fabs(Dqp) < fabs(Dqc))
      {
			  ans[var] = Dqp;
      }
		  else
      {
			  ans[var] = Dqc;
      }
	  }

  }

}
