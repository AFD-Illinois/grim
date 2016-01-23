#include "../problem.h"
#include "torus.h"

#if (CONDUCTION || VISCOSITY)
  void setDiffusionCoefficients(const struct geometry geom[ARRAY_ARGS 1],
                                struct fluidElement elem[ARRAY_ARGS 1]
                               )
  {
    #if (VISCOSITY)
      elem->nu  = 0.1;
    #endif
    #if (CONDUCTION)
      elem->chi = 0.1; 
    #endif
  }
#endif

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  VecSetRandom(ts->primN.vec, NULL);

//  startFillingVecGhost(&ts->primN);
//
//  finishFillingVecGhost(&ts->primN);
  setPointerToVec(&ts->primN);
  setPointerToVec(&ts->connection);
  setPointerToVec(&ts->conservedVarsN);
  setPointerToVec(&ts->sources);

  REAL ***ptr __attribute__((aligned(64)));
  DMDAVecGetArrayDOF(ts->primN.dm, ts->primN.vec, &ptr);

  struct geometry geom;
  REAL XCoords[NDIM] = {0., 10., .5, 3.};
  setGeometry(XCoords, &geom);
//  if ( ((unsigned long)ptr & 63) == 0)
//  {
//    PetscPrintf(PETSC_COMM_WORLD, "Ptr is aligned to 64 byte boundary\n");
//  } else
//  {
//    PetscPrintf(PETSC_COMM_WORLD, "Ptr is NOT aligned to 64 byte boundary\n");
//  }

  for (int n=0; n<250000; n++)
  {
    LOOP_OVER_TILES(&ts->primN)
    {
      for (int jInTile=0; jInTile<TILE_SIZE_X2; jInTile++)
      {
        REAL primVars[DOF][TILE_SIZE_X1]          __attribute__((aligned(64)));
        REAL primVarsTranspose[TILE_SIZE_X1][DOF] __attribute__((aligned(64)));
        REAL consVarsTranspose[TILE_SIZE_X1][DOF] __attribute__((aligned(64)));
        REAL gamma[TILE_SIZE_X1]                  __attribute__((aligned(64)));
        REAL uCon[NDIM][TILE_SIZE_X1]             __attribute__((aligned(64)));
        REAL uCov[NDIM][TILE_SIZE_X1]             __attribute__((aligned(64)));
        REAL bCon[NDIM][TILE_SIZE_X1]             __attribute__((aligned(64)));
        REAL bCov[NDIM][TILE_SIZE_X1]             __attribute__((aligned(64)));
        REAL pressure[TILE_SIZE_X1]               __attribute__((aligned(64)));
        REAL bSqr[TILE_SIZE_X1]                   __attribute__((aligned(64)));
        REAL moments[20][TILE_SIZE_X1]            __attribute__((aligned(64)));

        const int j = jInTile + jTile*TILE_SIZE_X2;

        #pragma vector aligned
        for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
        {
          const int i = iInTile + iTile*TILE_SIZE_X1;
          
          #pragma vector aligned
          for (int var=0; var<DOF; var++)
          {
            primVars[var][iInTile] = ts->primN.ptr[j][i][var];
          }
        }

        #pragma vector aligned
        for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
        {
          gamma[iInTile] = 
            sqrt(1 + geom.gCov[1][1]*primVars[U1][iInTile]*primVars[U1][iInTile]
                   + geom.gCov[2][2]*primVars[U2][iInTile]*primVars[U2][iInTile]
                   + geom.gCov[3][3]*primVars[U3][iInTile]*primVars[U3][iInTile]

                   + 2*(  geom.gCov[1][2]*primVars[U1][iInTile]*primVars[U2][iInTile]
                        + geom.gCov[1][3]*primVars[U1][iInTile]*primVars[U3][iInTile] 
                        + geom.gCov[2][3]*primVars[U2][iInTile]*primVars[U3][iInTile]
                       )
                );

          uCon[0][iInTile] = gamma[iInTile]/geom.alpha;

          uCon[1][iInTile] =  primVars[U1][iInTile]
                            - gamma[iInTile]*geom.gCon[0][1]*geom.alpha;

          uCon[2][iInTile] =  primVars[U2][iInTile]
                            - gamma[iInTile]*geom.gCon[0][2]*geom.alpha;

          uCon[3][iInTile] =  primVars[U3][iInTile]
                            - gamma[iInTile]*geom.gCon[0][3]*geom.alpha;
        }

        for (int mu=0; mu<NDIM; mu++)
        {
          for (int nu=0; nu<NDIM; nu++)
          {
            #pragma vector aligned
            #pragma simd
            for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
            { 
              uCov[mu][iInTile] += geom.gCov[mu][nu] * uCon[nu][iInTile];
            }
          }
        }

        #pragma vector aligned
        for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
        {
          bCon[0][iInTile] =  primVars[B1][iInTile]*uCov[1][iInTile]
                            + primVars[B2][iInTile]*uCov[2][iInTile]
                            + primVars[B3][iInTile]*uCov[3][iInTile];
  
          bCon[1][iInTile] =
              (primVars[B1][iInTile] + bCon[0][iInTile]*uCon[1][iInTile]
              )/uCon[0][iInTile];

          bCon[2][iInTile] = 
              (primVars[B2][iInTile] + bCon[0][iInTile]*uCon[2][iInTile]
              )/uCon[0][iInTile];

          bCon[3][iInTile] =
              (primVars[B3][iInTile] + bCon[0][iInTile]*uCon[3][iInTile]
              )/uCon[0][iInTile];
        }

        for (int mu=0; mu<NDIM; mu++)
        {
          for (int nu=0; nu<NDIM; nu++)
          {
            #pragma vector aligned
            #pragma simd
            for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
            { 
              bCov[mu][iInTile] += geom.gCov[mu][nu] * bCon[nu][iInTile];
            }
          }
        }

//        #pragma vector aligned
//        for (int iInTile=0; iInTile<TILE_SIZE_X1; iInTile++)
//        { 
//          pressure[iInTile] = (ADIABATIC_INDEX-1.)*primVars[UU][iInTile];
//          bSqr[iInTile]     =   bCon[0][iInTile]*bCov[0][iInTile]
//                              + bCon[1][iInTile]*bCov[1][iInTile]
//                              + bCon[2][iInTile]*bCov[2][iInTile]
//                              + bCon[3][iInTile]*bCov[3][iInTile];
//        }
        
        for (int mu=0; mu<1; mu++)
        {
          #pragma vector aligned
          for (int iInTile=0; iInTile < TILE_SIZE_X1; iInTile++)
          {
            moments[N_UP(mu)][iInTile] = 
              primVars[RHO][iInTile]*uCon[mu][iInTile];
          }

          for (int nu=0; nu<NDIM; nu++)
          {
            #pragma vector aligned
            for (int iInTile=0; iInTile < TILE_SIZE_X1; iInTile++)
            {
              REAL pressureScalar = (ADIABATIC_INDEX-1.)*primVars[UU][iInTile];
              REAL bSqrScalar     =   bCon[0][iInTile]*bCov[0][iInTile]
                              + bCon[1][iInTile]*bCov[1][iInTile]
                              + bCon[2][iInTile]*bCov[2][iInTile]
                              + bCon[3][iInTile]*bCov[3][iInTile];

              moments[T_UP_DOWN(mu,nu)][iInTile] =   
                            (  primVars[RHO][iInTile] + primVars[UU][iInTile]
                             + pressureScalar + bSqrScalar
                            ) * uCon[mu][iInTile] * uCov[nu][iInTile]

                          + (pressure[iInTile] + 0.5*bSqr[iInTile])*DELTA(mu, nu)

                          - bCon[mu][iInTile] * bCov[nu][iInTile]
                          ;
            }
          }
        }

        REAL g = sqrt(-geom.gDet);

        #pragma vector aligned
        #pragma ivdep
        for (int iInTile=0; iInTile < TILE_SIZE_X1; iInTile++)
        {
          const int i = iInTile + iTile*TILE_SIZE_X1;

          ts->conservedVarsN.ptr[j][i][RHO] = 
            g*moments[N_UP(0)][iInTile];

          ts->conservedVarsN.ptr[j][i][UU]  =  
            g*moments[T_UP_DOWN(0, 0)][iInTile] + g*moments[N_UP(0)][iInTile];

          ts->conservedVarsN.ptr[j][i][U1]  =
            g*moments[T_UP_DOWN(0, 1)][iInTile];

          ts->conservedVarsN.ptr[j][i][U2]  = 
            g*moments[T_UP_DOWN(0, 2)][iInTile];

          ts->conservedVarsN.ptr[j][i][U3]  = 
            g*moments[T_UP_DOWN(0, 3)][iInTile];

          ts->conservedVarsN.ptr[j][i][B1] = 
            g*(  bCon[1][iInTile]*uCon[0][iInTile] 
               - bCon[0][iInTile]*uCon[1][iInTile]
              );

          ts->conservedVarsN.ptr[j][i][B2] = 
            g*(  bCon[2][iInTile]*uCon[0][iInTile] 
               - bCon[0][iInTile]*uCon[2][iInTile]
              );

          ts->conservedVarsN.ptr[j][i][B3] = 
            g*(  bCon[3][iInTile]*uCon[0][iInTile] 
               - bCon[0][iInTile]*uCon[3][iInTile]
              );
        }

      }

    }

  }

  restorePointerToVec(&ts->primN);
  restorePointerToVec(&ts->connection);
  restorePointerToVec(&ts->conservedVarsN);
  restorePointerToVec(&ts->sources);
}

//void applyFloor(const int iTile, const int jTile,
//                const int X1Start, const int X2Start,
//                const int X1Size, const int X2Size,
//                REAL primTile[ARRAY_ARGS TILE_SIZE])
//{
//}
//
//void applyAdditionalProblemSpecificBCs
//(
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  const struct problemData problemSpecificData[ARRAY_ARGS 1],
//  REAL primTile[ARRAY_ARGS TILE_SIZE]
//)
//{
//}
//
//void applyProblemSpecificFluxFilter
//(
//  const int iTile, const int jTile,
//  const int X1Start, const int X2Start,
//  const int X1Size, const int X2Size,
//  const struct problemData problemSpecificData[ARRAY_ARGS 1],
//  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
//  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
//)
//{
//}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
}
