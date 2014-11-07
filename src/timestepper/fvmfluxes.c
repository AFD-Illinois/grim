#include "timestepper.h"

void computeFluxesOverTile(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                           const int iTile, const int jTile,
                           const int X1Start, const int X2Start,
                           const int X1Size, const int X2Size,
                           REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
                           REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE])
{
  REAL primVarsLeft[TILE_SIZE], primVarsRight[TILE_SIZE];
  REAL fluxTileLeft[TILE_SIZE], fluxTileRight[TILE_SIZE];
  REAL conservedVarsTileLeft[TILE_SIZE], conservedVarsTileRight[TILE_SIZE];

  reconstruct(primTile, X1,
              iTile, jTile, 
              X1Start, X2Start,
              X1Size, X2Size,
              primVarsLeft, primVarsRight);

  /* Requires NG>=3 */
  LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    REAL XCoords[NDIM];
    getXCoords(&zone, FACE_X1, XCoords);
    struct geometry geom; setGeometry(XCoords, &geom);
    struct fluidElement elem;
    setFluidElement(&primVarsLeft[INDEX_TILE(&zone, 0)], &geom, &elem);
  
    computeFluxes(&elem, &geom, 1,
                  &fluxTileLeft[INDEX_TILE(&zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileLeft[INDEX_TILE(&zone, 0)]);

    getXCoords(&zone, FACE_X1_PLUS_ONE, XCoords);
    setGeometry(XCoords, &geom);
    setFluidElement(&primVarsRight[INDEX_TILE(&zone, 0)], &geom, &elem);
  
    computeFluxes(&elem, &geom, 1,
                  &fluxTileRight[INDEX_TILE(&zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileRight[INDEX_TILE(&zone, 0)]);
  }
        
  LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    REAL XCoords[NDIM];
    getXCoords(&zone, FACE_X1, XCoords);
    struct geometry geom; setGeometry(XCoords, &geom);

    riemannSolver(&fluxTileRight[INDEX_TILE_MINUS_ONE_X1(&zone, 0)],
                  &fluxTileLeft[INDEX_TILE(&zone, 0)],
                  &conservedVarsTileRight[INDEX_TILE_MINUS_ONE_X1(&zone, 0)],
                  &conservedVarsTileLeft[INDEX_TILE(&zone, 0)],
                  &primVarsRight[INDEX_TILE_MINUS_ONE_X1(&zone, 0)],
                  &primVarsLeft[INDEX_TILE(&zone, 0)],
                  &geom, 1, &fluxX1Tile[INDEX_TILE(&zone, 0)]);
  }

#if (COMPUTE_DIM==2)
  reconstruct(primTile, X2,
              iTile, jTile, 
              X1Start, X2Start,
              X1Size, X2Size,
              primVarsLeft, primVarsRight);

  LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    REAL XCoords[NDIM];
    getXCoords(&zone, FACE_X2, XCoords);
    struct geometry geom; setGeometry(XCoords, &geom);
    struct fluidElement elem;
    setFluidElement(&primVarsLeft[INDEX_TILE(&zone, 0)], &geom, &elem);
  
    computeFluxes(&elem, &geom, 2,
                  &fluxTileLeft[INDEX_TILE(&zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileLeft[INDEX_TILE(&zone, 0)]);

    getXCoords(&zone, FACE_X2_PLUS_ONE, XCoords);
    setGeometry(XCoords, &geom);
    setFluidElement(&primVarsRight[INDEX_TILE(&zone, 0)], &geom, &elem);
  
    computeFluxes(&elem, &geom, 2,
                  &fluxTileRight[INDEX_TILE(&zone, 0)]);
    computeFluxes(&elem, &geom, 0,
                  &conservedVarsTileRight[INDEX_TILE(&zone, 0)]);
  }
        
  LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    REAL XCoords[NDIM];
    getXCoords(&zone, FACE_X2, XCoords);
    struct geometry geom; setGeometry(XCoords, &geom);

    riemannSolver(&fluxTileRight[INDEX_TILE_MINUS_ONE_X2(&zone, 0)],
                  &fluxTileLeft[INDEX_TILE(&zone, 0)],
                  &conservedVarsTileRight[INDEX_TILE_MINUS_ONE_X2(&zone, 0)],
                  &conservedVarsTileLeft[INDEX_TILE(&zone, 0)],
                  &primVarsRight[INDEX_TILE_MINUS_ONE_X2(&zone, 0)],
                  &primVarsLeft[INDEX_TILE(&zone, 0)],
                  &geom, 2, &fluxX2Tile[INDEX_TILE(&zone, 0)]);
  }

  /* Flux CT */
  REAL emf[TILE_SIZE];
  LOOP_INSIDE_TILE(0, TILE_SIZE_X1+2, 0, TILE_SIZE_X2+2)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    emf[INDEX_TILE(&zone, 0)] = 
                  0.25*(  fluxX1Tile[INDEX_TILE(&zone, B2)]
                        + fluxX1Tile[INDEX_TILE_MINUS_ONE_X2(&zone, B2)]
                        - fluxX2Tile[INDEX_TILE(&zone, B1)]
                        - fluxX2Tile[INDEX_TILE_MINUS_ONE_X1(&zone, B1)]
                       );
          
  }

  LOOP_INSIDE_TILE(0, TILE_SIZE_X1+1, 0, TILE_SIZE_X2+1)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);

    fluxX1Tile[INDEX_TILE(&zone, B1)] = 0.;

    fluxX1Tile[INDEX_TILE(&zone, B2)] = 
                                  0.5*(  emf[INDEX_TILE(&zone, 0)]
                                       + emf[INDEX_TILE_PLUS_ONE_X2(&zone, 0)]
                                      );

    fluxX2Tile[INDEX_TILE(&zone, B1)] = 
                                 -0.5*(  emf[INDEX_TILE(&zone, 0)]
                                       + emf[INDEX_TILE_PLUS_ONE_X1(&zone, 0)]
                                      );

    fluxX2Tile[INDEX_TILE(&zone, B2)] = 0.;

  }                                         
#endif /* COMPUTE_DIM==2 */

}
