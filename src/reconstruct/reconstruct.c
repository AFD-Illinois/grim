#include "reconstruct.h"

/* Taken from HARM */
void slopeLim(const REAL left[ARRAY_ARGS DOF],
              const REAL mid[ARRAY_ARGS DOF],
              const REAL right[ARRAY_ARGS DOF],
              REAL ans[ARRAY_ARGS DOF]
             )
{
  for (int var=0; var<DOF; var++)
  {
    #if (RECONSTRUCTION == MONOTONIZED_CENTRAL)
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
    #elif (RECONSTRUCTION == MIN_MOD)
		  REAL Dqm = (mid[var] - left[var]);
  		REAL Dqp = (right[var] - mid[var]);
	  	REAL s = Dqm*Dqp;
		  if(s <= 0.) ans[var] = 0.;
  		else if(fabs(Dqm) < fabs(Dqp)) ans[var] = Dqm;
	  	else ans[var] = Dqp;
    #endif
  }
}

void reconstruct(const REAL primTile[ARRAY_ARGS TILE_SIZE(DOF)],
                 const struct gridTile tile[ARRAY_ARGS 1],
                 const int dir,
                 REAL primTileLeft[ARRAY_ARGS TILE_SIZE(DOF)],
                 REAL primTileRight[ARRAY_ARGS TILE_SIZE(DOF)]
                )
{
  int indexMinus2, indexMinus1, index0, indexPlus1, indexPlus2;

  if (dir==1)
  {
    indexMinus2 = INDEX_TILE_OFFSET(-2, 0, 0, zone, 0);
    indexMinus1 = INDEX_TILE_OFFSET(-1, 0, 0, zone, 0);
    index0      = INDEX_TILE_OFFSET( 0, 0, 0, zone, 0);
    indexPlus1  = INDEX_TILE_OFFSET( 1, 0, 0, zone, 0);
    indexPlus2  = INDEX_TILE_OFFSET( 2, 0, 0, zone, 0);
  }
  else if (dir==2)
  {
    indexMinus2 = INDEX_TILE_OFFSET(0, -2, 0, zone, 0);
    indexMinus1 = INDEX_TILE_OFFSET(0, -1, 0, zone, 0);
    index0      = INDEX_TILE_OFFSET(0,  0, 0, zone, 0);
    indexPlus1  = INDEX_TILE_OFFSET(0,  1, 0, zone, 0);
    indexPlus2  = INDEX_TILE_OFFSET(0,  2, 0, zone, 0);
  }
  else if (dir==3)
  {
    indexMinus2 = INDEX_TILE_OFFSET(0, 0, -2, zone, 0);
    indexMinus1 = INDEX_TILE_OFFSET(0, 0, -1, zone, 0);
    index0      = INDEX_TILE_OFFSET(0, 0,  0, zone, 0);
    indexPlus1  = INDEX_TILE_OFFSET(0, 0,  1, zone, 0);
    indexPlus2  = INDEX_TILE_OFFSET(0, 0,  2, zone, 0);
  }

  LOOP_INSIDE_TILE(-1, TILE_SIZE_X1+1,
                   -1, TILE_SIZE_X2+1,
                   -1, TILE_SIZE_X3+1
                  )
  {
    struct gridZone zone;
    setGridZone(iInTile, jInTile, kInTile, tile, &zone);

    #if (RECONSTRUCTION == MONOTONIZED_CENTRAL || \
         RECONSTRUCTION == MIN_MOD \
        )
      REAL slope[DOF];

      slopeLim(&primTile[indexMinus1], 
               &primTile[index0], 
               &primTile[indexPlus1],
               slope
              );

      for (int var=0; var<DOF; var++)
      {
        /* Left Edge */
        primTileLeft[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE(&zone, var)] - 0.5*slope[var];

        /* Right Edge */
        primTileRight[INDEX_TILE(&zone, var)] =
          primTile[INDEX_TILE(&zone, var)] + 0.5*slope[var];
      }
    #endif
  }
}
