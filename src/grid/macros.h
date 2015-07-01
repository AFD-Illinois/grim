#ifndef GRIM_GRID_MACROS_H_
#define GRIM_GRID_MACROS_H_

#define CORNER           (0)
#define CENTER           (1)
#define FACE_X1          (2)
#define FACE_X1_PLUS_ONE (3)
#define FACE_X2          (4)
#define FACE_X2_PLUS_ONE (5)
#define FACE_X3          (6)
#define FACE_X3_PLUS_ONE (7)

#define TILE_SIZE_X3     (1)

/* Entire domain is divided into three parts:
 * 1) Global -- The full domain without ghost zones. Reading the global domain
 *              will give the part of the domain that the current processor can
 *              read. Using this means that the associated array does not need
 *              any communication. Used for arrays associated with global Petsc
 *              Vecs.
 *
 * 2) Local -- The part of the full domain that the current processor can read
 *             with ghost zones. Arrays marked local need communication. Used
 *             for arrays associated with local Petsc Vecs.
 *
 * 3) Tile -- Small 2D chunk of the local/global array that can fit inside the
 *            cache. Used for arrays created on our own. Each tile is described
 *            by (iTile, jTile) and a k index denoting its location in the X3
 *            direction. In 3D, we stream through the X3 direction, while tile
 *            blocking in the 2D X1-X2 plane.
 *
 * The prototypes for the following macros are the same for COMPUTE_DIM=1,2,3.
 * This allows the code to be dimension agnostic.
 *
 * Loop over the local grid is as follows:
 *
 * struct gridData *grid;
 * LOOP_OVER_TILES(grid)
 * {
 *   REAL tile[TILE_SIZE(DOF)];
 *   LOOP_INSIDE_TILE(0, DOF, -NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
 *   {
 *     struct gridZone zone;
 *     setZone(iInTile, jInTile, 
 *             iTile,   jTile,
 *             kGlobal, &zone
 *            );
 *
 *     Do computation over the tile
 *               .
 *               .
 *               .
 *   }
 *
 * }
 *
 * X1Size, X2Size, and X3 Size are determined using DMDAGetCorners() in
 * timeStepper
 *
 * */
#if (COMPUTE_DIM==1)

  /* Total size needed for a tile with {\tt numVar} variables */
  #define TILE_SIZE(numVar) ((TILE_SIZE_X1+2*NG)*(numVar))

  /* Returns the index of the zone */
  #define INDEX_TILE(zone,var) \
    ((zone)->iInTile + NG + (TILE_SIZE_X1+2*NG)*(var))

  /* Returns the index of the zone with an offset of (iOffset, jOffset, kOffset) */
  #define INDEX_TILE_OFFSET(iOffset,jOffset,kOffset,zone,var) \
    ((zone)->iInTile + (iOffset) + NG + (TILE_SIZE_X1+2*NG)*(var))

  /* Macro needed to access Petsc's global Vecs */
  #define INDEX_GRID(grid,zone,var) ((grid)->ptr[(zone)->iGlobal][var])

  /* Macro needed to access Petsc's local Vecs. */
  #define INDEX_GRID_GHOST(grid,zone,var) \
    ((grid)->ptrGhost[(zone)->iGlobal][var])

  /* Macro to define pointers that are needed to access Petsc Vecs */
  #define POINTER_TO_VEC(ptr) REAL **ptr /* [N1][N_VAR] */

#elif (COMPUTE_DIM==2)

  /* Total size needed for a tile with {\tt numVar} variables */
  #define TILE_SIZE(numVar) ((TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*(numVar))

  /* Returns the index of the zone */
  #define INDEX_TILE(zone,var) \
    ( (zone)->iInTile + NG + (TILE_SIZE_X1+2*NG)\
      * ( \
          (zone)->jInTile + NG + (TILE_SIZE_X2+2*NG)*(var)\
        ) \
    )

  /* Returns the index of the zone with an offset of (iOffset, jOffset, kOffset) */
  #define INDEX_TILE_OFFSET(iOffset,jOffset,kOffset,zone,var) \
    ( (zone)->iInTile + (iOffset) + NG + (TILE_SIZE_X1+2*NG)\
      * ( \
          (zone)->jInTile + (jOffset) + NG + (TILE_SIZE_X2+2*NG)*(var)\
        ) \
    )

  /* Macro needed to access Petsc's global Vecs */
  #define INDEX_GRID(grid,zone,var) \
    ((grid)->ptr[(zone)->jGlobal][(zone)->iGlobal][var])

  /* Macro needed to access Petsc's local Vecs */
  #define INDEX_GRID_GHOST(grid,zone,var) \
    ((grid)->ptrGhost[(zone)->jGlobal][(zone)->iGlobal][var])

  /* Macro to define pointers that are needed to access Petsc Vecs */
  #define POINTER_TO_VEC(ptr) REAL ***ptr /* [N2][N1][N_VAR] */

#elif (COMPUTE_DIM==3)
  /* We use the tile blocking technique from 
   * "3.5 D Blocking Optimization for Stencil Computations on Modern CPUs and
   * GPUs" -- A. Nguyen et. al., 2010 */

  /* Total size needed for a tile with {\tt numVar} variables. In 3D we allocate
   * (1+2*NG) elements in the X3 direction.*/
  #define TILE_SIZE(numVar) \
    ((TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*(1+2*NG)*(numVar))

  /* Returns the index of the zone */
  #define INDEX_TILE(zone,var) \
    ( (zone)->iInTile + NG + (TILE_SIZE_X1+2*NG)\
      * ( \
          (zone)->jInTile + NG + (TILE_SIZE_X2+2*NG)*(NG + (1+2*NG)*(var) )\
        ) \
    )

  /* Returns the index of the zone with an offset of (iOffset, jOffset, kOffset) */
  #define INDEX_TILE_OFFSET(iOffset,jOffset,kOffset,zone,var) \
    ( (zone)->iInTile + (iOffset) + NG + (TILE_SIZE_X1+2*NG)\
      * ( \
          (zone)->jInTile + (jOffset) + NG + (TILE_SIZE_X2+2*NG) \
         *((kOffset) + NG + (1+2*NG)*(var) )\
        ) \
    )

  /* Macro needed to access Petsc's global and local Vecs */
  #define INDEX_PETSC(ptr,zone,var) \
    (ptr[(zone)->kGlobal][(zone)->jGlobal][(zone)->iGlobal][var])

  /* Macro to define pointers that are needed to access Petsc Vecs */
  #define POINTER_TO_VEC(ptr) REAL ****ptr /* [N3][N2][N1][N_VAR] */

#endif

/* Macro to loop over tiles */
#if (USE_OPENMP)

  #define LOOP_OVER_TILES(grid) \
    _Pragma("omp parallel for") \
    for (int kGlobal = (grid)->kLocalStart; \
             kGlobal < (grid)->kLocalStart + (grid)->kLocalSize; kGlobal++) \
      for (int jTile = 0; jTile < (grid)->numTilesX2; jTile++) \
        for (int iTile = 0; iTile < (grid)->numTilesX1; iTile++)

#else

  #define LOOP_OVER_TILES(grid) \
    for (int kGlobal = (grid)->kLocalStart; \
             kGlobal < (grid)->kLocalStart + (grid)->kLocalSize; kGlobal++) \
      for (int jTile = 0; jTile < (grid)->numTilesX2; jTile++) \
        for (int iTile = 0; iTile < (grid)->numTilesX1; iTile++)

#endif

/* Macro to loop inside a tile */
#if (COMPUTE_DIM==1 || COMPUTE_DIM==2)
  #define LOOP_INSIDE_TILE(iInTileStart,  iInTileEnd, \
                           jInTileStart,  jInTileEnd, \
                           kInTileStart,  kInTileEnd \
                          ) \
    for (int kInTile = 0; kInTile < 1; kInTile++ ) \
      for (int jInTile = (jInTileStart); jInTile < (jInTileEnd); jInTile++ ) \
        _Pragma("vector aligned nontemporal") \
        _Pragma("simd") \
        for (int iInTile = (iInTileStart); iInTile < (iInTileEnd); iInTile++ )

#elif (COMPUTE_DIM==3)

  #define LOOP_INSIDE_TILE(iInTileStart,  iInTileEnd, \
                           jInTileStart,  jInTileEnd, \
                           kInTileStart,  kInTileEnd \
                          ) \
    for (int kInTile = (kInTileStart); kInTile < (kInTileEnd); kInTile++ ) \
      for (int jInTile = (jInTileStart); jInTile < (jInTileEnd); jInTile++ ) \
        _Pragma("vector aligned nontemporal") \
        _Pragma("simd") \
        for (int iInTile = (iInTileStart); iInTile < (iInTileEnd); iInTile++ )
#endif

#endif /* GRIM_GRID_MACROS_H_ */
