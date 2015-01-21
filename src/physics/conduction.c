#include "physics.h"

void addConductionSourceTermsToResidual
(
  const REAL primTile[ARRAY_ARGS TILE_SIZE],
  ARRAY(primGlobal), ARRAY(primHalfStepGlobal), ARRAY(primOldGlobal),
  ARRAY(connectionGlobal),
  ARRAY(gradTGlobal), ARRAY(graduConGlobal), 
  ARRAY(graduConHigherOrderTerm1Global),
  REAL dt,
  int computeOldSourceTermsAndOldDivOfFluxes,
  int computeDivOfFluxAtTimeN,
  int computeDivOfFluxAtTimeNPlusHalf,
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  ARRAY(residualGlobal)
)
{
#if (CONDUCTION)
  LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
  {
    struct gridZone zoneCenter;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zoneCenter);
    
    /* Values needed to compute spatial derivative.  Depending on where the
     * spatial derivatives are to be computed (t=n or t=n+1/2), the values of
     * primTile is set (already done, outside this routine) */

    if (computeOldSourceTermsAndOldDivOfFluxes)
    {
      REAL gradT[COMPUTE_DIM];
      REAL graduCon[COMPUTE_DIM*NDIM];
      REAL graduConHigherOrderTerm1[COMPUTE_DIM];
      computeConductionSpatialGradientTerms
        (primTile, 
         iTile, jTile, 
         iInTile, jInTile,
         X1Start, X2Start,
         X1Size, X2Size, 
         gradT, graduCon, graduConHigherOrderTerm1
        );
      
      INDEX_PETSC(gradTGlobal, &zoneCenter, 0) = gradT[0];

      for (int mu=0; mu<NDIM; mu++)
      {
        INDEX_PETSC(graduConGlobal, &zoneCenter, mu) = graduCon[mu];
      }

      INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 0) = 
        graduConHigherOrderTerm1[0];

      #if (COMPUTE_DIM==2)
        INDEX_PETSC(gradTGlobal, &zoneCenter, 1) = gradT[1];
        
        for (int mu=0; mu<NDIM; mu++)
        {
          INDEX_PETSC(graduConGlobal, &zoneCenter, mu+NDIM) = graduCon[mu+NDIM];
        }

        INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 1) = 
          graduConHigherOrderTerm1[1];
      #endif
    } 
    else
    {
      #if (TIME_STEPPING==IMPLICIT)
        REAL gradT[COMPUTE_DIM];
        REAL graduCon[COMPUTE_DIM*NDIM];
        REAL graduConHigherOrderTerm1[COMPUTE_DIM];
        computeConductionSpatialGradientTerms
          (primTile, 
           iTile, jTile, 
           iInTile, jInTile,
           X1Start, X2Start,
           X1Size, X2Size, 
           gradT, graduCon, graduConHigherOrderTerm1
          );
      #endif

      /* Values needed to compute time derivative */
      REAL T    =   (ADIABATIC_INDEX-1.)
                  * INDEX_PETSC(primGlobal, &zoneCenter, UU)
                  / INDEX_PETSC(primGlobal, &zoneCenter, RHO);

      REAL TOld =   (ADIABATIC_INDEX-1.)
                  * INDEX_PETSC(primOldGlobal, &zoneCenter, UU)
                  / INDEX_PETSC(primOldGlobal, &zoneCenter, RHO);

      REAL dT[NDIM];

      dT[0] = (T - TOld)/dt;

      #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
        dT[1] = INDEX_PETSC(gradTGlobal, &zoneCenter, 0);

        dT[2] = 0.; dT[3] = 0.;
    
        #if (COMPUTE_DIM==2)
          dT[2] = INDEX_PETSC(gradTGlobal, &zoneCenter, 1);
        #endif

      #elif (TIME_STEPPING == IMPLICIT)
        dT[1] = 0.5*(  INDEX_PETSC(gradTGlobal, &zoneCenter, 0)
                     + gradT[0]
                    );

        dT[2] = 0.; dT[3] = 0.;
    
        #if (COMPUTE_DIM==2)
          dT[2] = 0.5*(  INDEX_PETSC(gradTGlobal, &zoneCenter, 1)
                       + gradT[1]
                      );
        #endif

      #endif

      REAL aCon[NDIM];

      struct fluidElement elem, elemCenter, elemOld;
      struct geometry geomCenter;

      REAL XCoords[NDIM];
      getXCoords(&zoneCenter, CENTER, XCoords);
      setGeometry(XCoords, &geomCenter);
      setFluidElement(&INDEX_PETSC(primGlobal, &zoneCenter, 0), &geomCenter,
                      &elem);

      setFluidElement(&INDEX_PETSC(primOldGlobal, &zoneCenter, 0), &geomCenter,
                      &elemOld);

      if (computeDivOfFluxAtTimeN)
      {
        setFluidElement(&INDEX_PETSC(primOldGlobal, &zoneCenter, 0), 
                        &geomCenter, &elemCenter);
      } 
      else if (computeDivOfFluxAtTimeNPlusHalf)
      {
        setFluidElement(&INDEX_PETSC(primHalfStepGlobal, &zoneCenter, 0),
                        &geomCenter, &elemCenter);
      }

      /* Higher order term 1 
       * phi*(g*u^\mu)_{,\mu} */

      REAL dHigherOrderTerm1[NDIM];
      
      dHigherOrderTerm1[0] = 
        elemCenter.primVars[PHI] * sqrt(-geomCenter.gDet) 
      * (elem.uCon[0] - elemOld.uCon[0])/dt;
      
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)
        dHigherOrderTerm1[1] = 
          INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 0);

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
            INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 1);
        #endif
      #elif (TIME_STEPPING == IMPLICIT)
        dHigherOrderTerm1[1] = 
          0.5*(  INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 0)
               + graduConHigherOrderTerm1[0]
              );

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
            0.5*(  INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 1)
                 + graduConHigherOrderTerm1[1]
                );
        #endif

      #endif

      REAL higherOrderTerm1 = 0.;
      for (int mu=0; mu<NDIM; mu++)
      {
        higherOrderTerm1 += dHigherOrderTerm1[mu];
      }

      for (int mu=0; mu<NDIM; mu++)
      {
        aCon[mu]  = 0.;

        aCon[mu] +=   elemCenter.uCon[0]
                    * (elem.uCon[mu] - elemOld.uCon[mu])/dt;
      
        #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)
          aCon[mu] +=   elemCenter.uCon[1]
                      * INDEX_PETSC(graduConGlobal, &zoneCenter, mu);

          #if (COMPUTE_DIM==2)
            aCon[mu] +=   elemCenter.uCon[2]
                        * INDEX_PETSC(graduConGlobal, &zoneCenter, mu+NDIM);
          #endif

        #elif (TIME_STEPPING == IMPLICIT)
          aCon[mu] += 
            0.5*(  (elemOld.uCon[1]*INDEX_PETSC(graduConGlobal, &zoneCenter, mu) )
                 + (elem.uCon[1]*graduCon[mu])
                );

          #if (COMPUTE_DIM==2)
            aCon[mu] += 
              0.5*(  (elemOld.uCon[2]*INDEX_PETSC(graduConGlobal, &zoneCenter, mu+NDIM) )
                   + (elem.uCon[2]*graduCon[mu])
                  );
          #endif

        #endif


        #if (TIME_STEPPING == IMEX || TIME_STEPPING == IMPLICIT)
          for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int beta=0; beta<NDIM; beta++)
            {
              aCon[mu] += 
              INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, beta))
              * (   (elem.uCon[alpha]    * elem.uCon[beta])
                  + (elemOld.uCon[alpha] * elemOld.uCon[beta])
                )/2.;
            }
          }
        #elif (TIME_STEPPING == EXPLICIT)
          for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int beta=0; beta<NDIM; beta++)
            {
              aCon[mu] +=  
              INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, beta))
              * elemCenter.uCon[alpha] * elemCenter.uCon[beta];
            }
          }
        #endif
      }

      REAL aCov[NDIM];
      for (int mu=0; mu<NDIM; mu++)
      {
        aCov[mu] = 0.;
        for (int nu=0; nu<NDIM; nu++)
        {
          aCov[mu] += geomCenter.gCov[mu][nu]*aCon[nu];
        }
      }

      REAL qConEckart[NDIM];
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)
        REAL TCenter =   (ADIABATIC_INDEX-1.)
                       * elemCenter.primVars[UU]
                       / elemCenter.primVars[RHO];
      #elif (TIME_STEPPING == IMPLICIT)
        REAL T =   (ADIABATIC_INDEX-1.)
                 * elem.primVars[UU]
                 / elem.primVars[RHO];

        REAL TOld =   (ADIABATIC_INDEX-1.)
                    * elemOld.primVars[UU]
                    / elemOld.primVars[RHO];

        REAL TCenter = (T + TOld)/2.;
      #endif

      for (int mu=0; mu<NDIM; mu++)
      {
        qConEckart[mu] = 0.;
        for (int nu=0; nu<NDIM; nu++)
        {
          qConEckart[mu] +=
            -elemCenter.kappa
            * (  elemCenter.uCon[mu]*elemCenter.uCon[nu] 
               + geomCenter.gCon[mu][nu]
              )
            * (dT[nu] + TCenter*aCov[nu]);
        }
      }

      REAL bDotq = 0., bSqr, bCov[NDIM];
      conToCov(elemCenter.bCon, &geomCenter, bCov);
      bSqr = getbSqr(&elemCenter, &geomCenter);

      for (int mu=0; mu<NDIM; mu++)
      {
        bDotq += bCov[mu]*qConEckart[mu]/sqrt(bSqr);
      }

      REAL g = sqrt(-geomCenter.gDet);
      #if (TIME_STEPPING == EXPLICIT)

        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) += 
          - higherOrderTerm1
          - g*(elemCenter.primVars[PHI] - bDotq)/elemCenter.tau;

      #elif (TIME_STEPPING == IMEX || TIME_STEPPING == IMPLICIT)
  
        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) += 
          - higherOrderTerm1
          - g*( 0.5*(elem.primVars[PHI] + elemOld.primVars[PHI]) 
               - bDotq
              )/elemCenter.tau;

      #endif

    }
  }
#endif /* CONDUCTION */
}

void computeConductionSpatialGradientTerms
(
  const REAL primTile[ARRAY_ARGS TILE_SIZE],
  const int iTile, const int jTile,
  const int iInTile, const int jInTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  REAL gradT[COMPUTE_DIM],
  REAL graduCon[COMPUTE_DIM*NDIM],
  REAL graduConHigherOrderTerm1[COMPUTE_DIM]
)
{
#if (CONDUCTION)
  REAL XCoords[NDIM];
  struct gridZone zone, zoneCenter;
  struct geometry geom;
  struct fluidElement elem, elemCenter;
  REAL uConLeft[NDIM], uConCenter[NDIM], uConRight[NDIM];
  REAL gLeft, gCenter, gRight;

  setGridZone(iTile, jTile,
              iInTile, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zoneCenter);

  REAL TCenter  =   (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE(&zoneCenter, UU)]
                  / primTile[INDEX_TILE(&zoneCenter, RHO)];

  REAL TLeftX1  =   (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE_MINUS_ONE_X1(&zoneCenter, UU)]
                  / primTile[INDEX_TILE_MINUS_ONE_X1(&zoneCenter, RHO)];

  REAL TRightX1 =  (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE_PLUS_ONE_X1(&zoneCenter, UU)]
                  / primTile[INDEX_TILE_PLUS_ONE_X1(&zoneCenter, RHO)];

  gradT[0] = 
      slopeLimitedDerivative(TLeftX1, TCenter, TRightX1)/zoneCenter.dX1;

  #if (COMPUTE_DIM==2)
    REAL TLeftX2  =   (ADIABATIC_INDEX-1.)
                    * primTile[INDEX_TILE_MINUS_ONE_X2(&zoneCenter, UU)]
                    / primTile[INDEX_TILE_MINUS_ONE_X2(&zoneCenter, RHO)];

    REAL TRightX2 =   (ADIABATIC_INDEX-1.)
                    * primTile[INDEX_TILE_PLUS_ONE_X2(&zoneCenter, UU)]
                    / primTile[INDEX_TILE_PLUS_ONE_X2(&zoneCenter, RHO)];

    gradT[1] = 
      slopeLimitedDerivative(TLeftX2, TCenter, TRightX2)/zoneCenter.dX2;
  #endif

  /* uConCenter */
  getXCoords(&zoneCenter, CENTER, XCoords);
  setGeometry(XCoords, &geom); gCenter = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zoneCenter, 0)], &geom, &elem);
  elemCenter = elem;
  for (int mu=0; mu<NDIM; mu++)
  {
    uConCenter[mu] = elem.uCon[mu];
  }

  /* uConLeft */
  setGridZone(iTile-1, jTile,
              iInTile, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zone);
  getXCoords(&zone, CENTER, XCoords);
  setGeometry(XCoords, &geom); gLeft = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConLeft[mu] = elem.uCon[mu];
  }

  /* uConRight */
  setGridZone(iTile+1, jTile,
              iInTile, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zone);
  getXCoords(&zone, CENTER, XCoords);
  setGeometry(XCoords, &geom); gRight = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConRight[mu] = elem.uCon[mu];
  }

  for (int mu=0; mu<NDIM; mu++)
  {
    graduCon[mu] = 
      slopeLimitedDerivative(uConLeft[mu],
                             uConCenter[mu],
                             uConRight[mu])/zoneCenter.dX1;
  }

  graduConHigherOrderTerm1[0] = 
    elemCenter.primVars[PHI]
  * slopeLimitedDerivative(gLeft*uConLeft[1],
                           gCenter*uConCenter[1],
                           gRight*uConRight[1])/zoneCenter.dX1;

  #if (COMPUTE_DIM==2)
    /* uConLeft */
    setGridZone(iTile, jTile-1,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    getXCoords(&zone, CENTER, XCoords);
    setGeometry(XCoords, &geom); gLeft = sqrt(-geom.gDet);
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
    for (int mu=0; mu<NDIM; mu++)
    {
      uConLeft[mu] = elem.uCon[mu];
    }

    /* uConRight */
    setGridZone(iTile, jTile+1,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    getXCoords(&zone, CENTER, XCoords);
    setGeometry(XCoords, &geom); gRight = sqrt(-geom.gDet);
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
    for (int mu=0; mu<NDIM; mu++)
    {
      uConRight[mu] = elem.uCon[mu];
    }

    for (int mu=0; mu<NDIM; mu++)
    {
      graduCon[mu+NDIM] = 
        slopeLimitedDerivative(uConLeft[mu],
                               uConCenter[mu],
                               uConRight[mu])/zoneCenter.dX2;
    }

    graduConHigherOrderTerm1[1] = 
      elemCenter.primVars[PHI]
    * slopeLimitedDerivative(gLeft*uConLeft[2],
                             gCenter*uConCenter[2],
                             gRight*uConRight[2])/zoneCenter.dX2;
  #endif

#endif /* CONDUCTION */
}
