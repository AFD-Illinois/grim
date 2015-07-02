#include "timestepper.h"

void fluxCT(const struct gridTile tile[ARRAY_ARGS 1],
            REAL fluxesTile[ARRAY_ARGS COMPUTE_DIM][TILE_SIZE(DOF)]
           )
{
  #if (COMPUTE_DIM==1)
    return;
  #endif

  REAL emfX3[TILE_SIZE(1)];
  memset(emfX3, 0., sizeof(REAL[TILE_SIZE(1)]) );

  #if (COMPUTE_DIM==3)
    REAL emfX1[TILE_SIZE(1)];
    REAL emfX2[TILE_SIZE(1)];
    memset(emfX1, 0., sizeof(REAL[TILE_SIZE(1)]) );
    memset(emfX2, 0., sizeof(REAL[TILE_SIZE(1)]) );
  #endif

  LOOP_INSIDE_TILE(0, TILE_SIZE_X1+1,
                   0, TILE_SIZE_X2+1,
                   0, TILE_SIZE_X3+1
                  )
  {
    struct gridZone zone;
    setZone(iInTile, jInTile, kInTile, tile, &zone);

    emfX3[INDEX_TILE(&zone, 0)] = 
      0.25*(  fluxesTile[X1][INDEX_TILE(&zone, B2)]
            + fluxesTile[X1][INDEX_TILE_OFFSET( 0, -1,  0, &zone, B2)]
            - fluxesTile[X2][INDEX_TILE(&zone, B1)]
            - fluxesTile[X2][INDEX_TILE_OFFSET(-1,  0,  0, &zone, B1)]
           );

    #if (COMPUTE_DIM==3)
      emfX1[INDEX_TILE(&zone, 0)] = 
        0.25*(  fluxesTile[X2][INDEX_TILE(&zone, B3)]
              + fluxesTile[X2][INDEX_TILE_OFFSET( 0,  0, -1, &zone, B3)]
              - fluxesTile[X3][INDEX_TILE(&zone, B2)]
              - fluxesTile[X3][INDEX_TILE_OFFSET( 0, -1,  0, &zone, B2)]
             );

      emfX2[INDEX_TILE(&zone, 0)] = 
        0.25*(  fluxesTile[X3][INDEX_TILE(&zone, B1)]
              + fluxesTile[X3][INDEX_TILE_OFFSET(-1,  0,  0, &zone, B1)]
              - fluxesTile[X1][INDEX_TILE(&zone, B3)]
              - fluxesTile[X1][INDEX_TILE_OFFSET( 0,  0, -1, &zone, B3)]
             );
    #endif
  }

  struct gridZone zoneAtTileEdge;
  setZone(0, 0, 0, tile, &zoneAtTileEdge);

  memset(&fluxesTile[X1][INDEX_TILE(&zoneAtTileEdge, B1)], 
         0., sizeof(REAL[TILE_SIZE(1)])
        );
  memset(&fluxesTile[X2][INDEX_TILE(&zoneAtTileEdge, B2)], 
         0., sizeof(REAL[TILE_SIZE(1)])
        );
  memset(&fluxesTile[X3][INDEX_TILE(&zoneAtTileEdge, B3)], 
         0., sizeof(REAL[TILE_SIZE(1)])
        );

  LOOP_INSIDE_TILE(0, TILE_SIZE_X1+1,
                   0, TILE_SIZE_X2+1,
                   0, TILE_SIZE_X3+1
                  )
  {
    struct gridZone zone;
    setZone(iInTile, jInTile, kInTile, tile, &zone);

    fluxesTile[X1][INDEX_TILE(&zone, B2)] =
      0.5*(  emfX3[INDEX_TILE(&zone, 0)]
           + emfX3[INDEX_TILE_OFFSET(0, 1, 0, &zone, 0)]
          );

    fluxesTile[X2][INDEX_TILE(&zone, B1)] =
      -0.5*(  emfX3[INDEX_TILE(&zone, 0)]
            + emfX3[INDEX_TILE_OFFSET(1, 0, 0, &zone, 0)]
           );

    #if (COMPUTE_DIM==3)

    fluxesTile[X1][INDEX_TILE(&zone, B3)] =
      -0.5*(  emfX2[INDEX_TILE(&zone, 0)]
            + emfX2[INDEX_TILE_OFFSET(0, 0, 1, &zone, 0)]
           );
    
    fluxesTile[X2][INDEX_TILE(&zone, B3)] =
      0.5*(  emfX1[INDEX_TILE(&zone, 0)]
           + emfX1[INDEX_TILE_OFFSET(0, 0, 1, &zone, 0)]
          );

    fluxesTile[X3][INDEX_TILE(&zone, B1)] =
      0.5*(  emfX2[INDEX_TILE(&zone, 0)]
           + emfX2[INDEX_TILE_OFFSET(1, 0, 0, &zone, 0)]
          );

    fluxesTile[X3][INDEX_TILE(&zone, B2)] =
      -0.5*(  emfX1[INDEX_TILE(&zone, 0)]
            + emfX1[INDEX_TILE_OFFSET(0, 1, 0, &zone, 0)]
           );

    #endif
  }

}
