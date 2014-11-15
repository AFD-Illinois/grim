#include "../problem.h"
#include "torus_2d.h"

/* Calculate the constant angular momentum per unit inertial mass (l = u_phi *
 * u^t) for a given black hole spin and a radius of the accretion disk.  Eqn 3.8
 * of Fishbone and Moncrief, 1976 */
REAL lFishboneMoncrief(REAL a, REAL r, REAL theta)
{
  return sqrt(M/pow(r, 3.)) \
        *(  pow(r, 4.) + r*r*a*a - 2.*M*r*a*a \
          - a*sqrt(M*r)*(r*r - a*a) \
         )/ \
         (r*r - 3*M*r + 2.*a*sqrt(M*r));
}

REAL lnOfhTerm1(REAL a,
                REAL r, REAL theta, 
                REAL l)
{
  REAL Delta = computeDelta(a, r, theta);
  REAL Sigma = computeSigma(a, r, theta);
  REAL A     = computeA(a, r, theta);

  return 0.5*log( (1. + sqrt(1. + (4.*l*l*Sigma*Sigma*Delta)/ \
                                  (A*sin(theta)*A*sin(theta))
                            )
                  ) / (Sigma*Delta/A)
                );
}

REAL lnOfhTerm2(REAL a,
                REAL r, REAL theta, 
                REAL l)
{
  REAL Delta = computeDelta(a, r, theta);
  REAL Sigma = computeSigma(a, r, theta);
  REAL A     = computeA(a, r, theta);

  return -0.5*sqrt(1. + (4.*l*l*Sigma*Sigma*Delta) /
                        (A*A*sin(theta)*sin(theta))
                  );

}

REAL lnOfhTerm3(REAL a,
                REAL r, REAL theta, 
                REAL l)
{
  REAL A     = computeA(a, r, theta);

  return -2*a*M*r*l/A;
}

REAL computeDelta(REAL a, REAL r, REAL theta)
{
  return r*r - 2*M*r + a*a;
}

REAL computeSigma(REAL a, REAL r, REAL theta)
{
  return r*r + a*a*cos(theta)*cos(theta);
}

REAL computeA(REAL a, REAL r, REAL theta)
{
  REAL Delta = computeDelta(a, r, theta);

  return pow(r*r + a*a, 2.) - Delta*a*a*sin(theta)*sin(theta);
}


REAL computeLnOfh(REAL a, REAL r, REAL theta)
{
  REAL l = lFishboneMoncrief(a, R_PRESSURE_MAX, M_PI/2.);

  REAL term1 = lnOfhTerm1(a, r, theta, l);
  REAL term2 = lnOfhTerm2(a, r, theta, l);
  REAL term3 = lnOfhTerm3(a, r, theta, l);

  REAL term1InnerEdge = lnOfhTerm1(a, R_INNER_EDGE, M_PI/2., l);
  REAL term2InnerEdge = lnOfhTerm2(a, R_INNER_EDGE, M_PI/2., l);
  REAL term3InnerEdge = lnOfhTerm3(a, R_INNER_EDGE, M_PI/2., l);

  return  term1 + term2 + term3 \
        - term1InnerEdge - term2InnerEdge - term3InnerEdge;

}

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL rhoMax=0., uMax=0.;
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM], xCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      XTox(XCoords, xCoords);

      /* The Fishbone-Moncrief solution is a stable hydrodynamic disk with a
       * constant angular momentum per unit inertial mass (l = u_phi * u^t). The
       * value of l is chosen by first choosing the radius at which the pressure
       * is maximum in the disk. Then the solution is characterized by two
       * free parameters:
       * 1) R_PRESSURE_MAX : The location of the pressure maximum of the disk.
       * 2) R_INNER_EDGE : The location of the inner edge of the disk.
       *
       * In addition, one has the choice of the adiabat C in the equation of
       * state: P = C rho^Gamma
       */

      REAL r = xCoords[1], theta = xCoords[2];

      REAL lnOfh;
      if (r >= R_INNER_EDGE)
      {
        lnOfh = computeLnOfh(BH_SPIN, r, theta);
      } else
      {
        lnOfh = 1.;
      }

      /* Region outside the torus */
      if (lnOfh < 0. || r < R_INNER_EDGE)
      {
        primTile[INDEX_TILE(&zone, RHO)] = 1e-15;
        primTile[INDEX_TILE(&zone, UU)] = 1e-15;
        primTile[INDEX_TILE(&zone, U1)] = 0.;
        primTile[INDEX_TILE(&zone, U2)] = 0.;
        primTile[INDEX_TILE(&zone, U3)] = 0.;
      } else
      {
        REAL h = exp(lnOfh);

        /* Solve for rho using the definition of h = (rho + u + P)/rho where rho
         * here is the rest mass energy density and P = C * rho^Gamma */

        REAL rho = pow((h-1)*(ADIABATIC_INDEX-1.)/(ADIABAT*ADIABATIC_INDEX), 
                       1./(ADIABATIC_INDEX-1.)
                      );
        REAL u =  ADIABAT * pow(rho, ADIABATIC_INDEX)
                / (ADIABATIC_INDEX-1.);


        /* Fishbone-Moncrief u_phi is given in the Boyer-Lindquist coordinates.
         * Need to transform to (modified) Kerr-Schild */
        REAL A = computeA(BH_SPIN, r, theta);
        REAL Sigma = computeSigma(BH_SPIN, r, theta);
        REAL Delta = computeDelta(BH_SPIN, r, theta);
        REAL l = lFishboneMoncrief(BH_SPIN, r, theta);

			  REAL expOfMinus2Chi = Sigma*Sigma*Delta/(A*A*sin(theta)*sin(theta)) ;
        REAL uCovPhiBL = sqrt((-1. + sqrt(1. + 4*l*l*expOfMinus2Chi)
                               )/2.
                              );
        REAL uConPhiBL =   2.*BH_SPIN*r*sqrt(1. + uCovPhiBL*uCovPhiBL)
                          /sqrt(A*Sigma*Delta)
                        + sqrt(Sigma/A)*uCovPhiBL/sin(theta);

        REAL uConBL[NDIM];
        uConBL[0] = 0.;
        uConBL[1] = 0.;
        uConBL[2] = uConPhiBL;

        REAL gCovBL[NDIM][NDIM], gConBL[NDIM][NDIM];
        REAL transformBLToMKS[NDIM][NDIM];
      
        for (int alpha=0; alpha<NDIM; alpha++)
        {
          for (int beta=0; beta<NDIM; beta++)
          {
            gCovBL[alpha][beta] = 0.;
            gConBL[alpha][beta] = 0.;
            transformBLToMKS[alpha][beta] = 0.;
          }
        }

        REAL mu = 1 + BH_SPIN*BH_SPIN*cos(theta)*cos(theta)/(r*r);

        gCovBL[0][0] = -(1. - 2./(r*mu));
        gCovBL[0][3] = -2.*BH_SPIN*sin(theta)*sin(theta)/(r*mu);
        gCovBL[3][0] = gCovBL[0][3];
        gCovBL[1][1] = mu*r*r/Delta;
        gCovBL[2][2] = r*r*mu;
        gCovBL[3][3] = r*r*sin(theta)*sin(theta)*\
                       (1. + BH_SPIN*BH_SPIN/(r*r) +
                        2.*BH_SPIN*BH_SPIN*sin(theta)*sin(theta)/\
                        (r*r*r*mu)
                       );

        gConBL[0][0] = -1. -2.*(1 + BH_SPIN*BH_SPIN/(r*r))/(Delta*mu/r);
        gConBL[0][3] = -2.*BH_SPIN/(r*Delta*mu);
        gConBL[3][0] = gConBL[0][3];
        gConBL[1][1] = Delta/(r*r*mu);
        gConBL[2][2] = 1./(r*r*mu);
        gConBL[3][3] = (1. - 2./(r*mu))/(sin(theta)*sin(theta)*Delta);

        transformBLToMKS[0][0] = 1.;
        transformBLToMKS[1][1] = 1.;
        transformBLToMKS[2][2] = 1.;
        transformBLToMKS[3][3] = 1.;
        transformBLToMKS[0][1] = 2.*r/Delta;
        transformBLToMKS[3][1] = BH_SPIN/Delta; 

        /* Need to get uConBL[0] using u^mu u_mu = -1 */
        REAL AA = gCovBL[0][0];
        REAL BB = 2.*(gCovBL[0][1]*uConBL[1] +
                      gCovBL[0][2]*uConBL[2] +
                      gCovBL[0][3]*uConBL[3]
                     );
        REAL CC = 1. + gCovBL[1][1]*uConBL[1]*uConBL[1] +
                       gCovBL[2][2]*uConBL[2]*uConBL[2] +
                       gCovBL[3][3]*uConBL[3]*uConBL[3] +
                   2.*(gCovBL[1][2]*uConBL[1]*uConBL[2] +
                       gCovBL[1][3]*uConBL[1]*uConBL[3] +
                       gCovBL[2][3]*uConBL[2]*uConBL[3]);
  
        REAL discriminent = BB*BB - 4.*AA*CC;
        uConBL[0] = -(BB + sqrt(discriminent))/(2.*AA);

        REAL uConKS[NDIM];

        for (int alpha=0; alpha<NDIM; alpha++)
        {
          uConKS[alpha] = 0.;

          for (int beta=0; beta<NDIM; beta++)
          {
            uConKS[alpha] += transformBLToMKS[alpha][beta]*uConBL[beta];
          }
        }

        /* Finally get the four-velocity in the X coordinates, which is modified
         * Kerr-Schild */
        REAL uConMKS[NDIM];
        REAL rFactor = r;
        REAL hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*XCoords[2]);
        uConMKS[0] = uConKS[0];
        uConMKS[1] = uConKS[1]/rFactor;
        uConMKS[2] = uConKS[2]/hFactor;
        uConMKS[3] = uConKS[3];


        if (rho > rhoMax) rhoMax = rho;
        if (u > uMax) uMax = u;

        primTile[INDEX_TILE(&zone, RHO)] = rho;
        primTile[INDEX_TILE(&zone, UU)] = u;

        struct geometry geom; 
        setGeometry(XCoords, &geom);
        primTile[INDEX_TILE(&zone, U1)] =   
          uConMKS[1] + pow(geom.alpha, 2.)*geom.gCon[0][1]*uConMKS[0];
        primTile[INDEX_TILE(&zone, U2)] =   
          uConMKS[2] + pow(geom.alpha, 2.)*geom.gCon[0][2]*uConMKS[0];
        primTile[INDEX_TILE(&zone, U3)] =   
          uConMKS[3] + pow(geom.alpha, 2.)*geom.gCon[0][3]*uConMKS[0];

        primTile[INDEX_TILE(&zone, B1)] = 0.;
        primTile[INDEX_TILE(&zone, B2)] = 0.;
        primTile[INDEX_TILE(&zone, B3)] = 0.;
      }
  
    }

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      primTile[INDEX_TILE(&zone, RHO)] = 
        primTile[INDEX_TILE(&zone, RHO)]/rhoMax;
      primTile[INDEX_TILE(&zone, UU)] =
        primTile[INDEX_TILE(&zone, UU)]/rhoMax;
    }

    applyFloor(iTile, jTile, 
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);


  }
  
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);
}


void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    REAL XCoords[NDIM], xCoords[NDIM];
    getXCoords(&zone, CENTER, XCoords);
    XTox(XCoords, xCoords);

    REAL r = xCoords[1];

    REAL rhoFloor = RHO_FLOOR*pow(r, RHO_FLOOR_FALLOFF);
    REAL uFloor = UU_FLOOR*pow(r, UU_FLOOR_FALLOFF);

    if (rhoFloor < RHO_FLOOR_MIN)
    {
      rhoFloor = RHO_FLOOR_MIN;
    }

    if (uFloor < UU_FLOOR_MIN)
    {
      uFloor = UU_FLOOR_MIN;
    }

    if (primTile[INDEX_TILE(&zone, RHO)] < rhoFloor)
    {
      primTile[INDEX_TILE(&zone, RHO)] = rhoFloor;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < uFloor)
    {
      primTile[INDEX_TILE(&zone, UU)] = uFloor;
    }

    /* Inflow check at the inner radial boundary */
    if (zone.i < 0)
    {
      inflowCheck(&zone, XCoords, primTile);
    }

    /* Inflow check at the outer radial boundary */
    if (zone.i > N1-1)
    {
      inflowCheck(&zone, XCoords, primTile);
    }

  }

}

void problemDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}

void inflowCheck(const struct gridZone zone[ARRAY_ARGS 1],
                 const REAL XCoords[NDIM],
                 REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  struct geometry geom;
  setGeometry(XCoords, &geom);

  struct fluidElement elem;
  setFluidElement(&primTile[INDEX_TILE(zone, 0)], &geom, &elem);

  int inflowInnerBoundary = (zone->i < 0 && elem.uCon[1] > 0);
  int inflowOuterBoundary = (zone->i >= N1 && elem.uCon[1] < 0);

  if (inflowInnerBoundary || inflowOuterBoundary)
  {
    /* Set new radial velocity to be zero */
    primTile[INDEX_TILE(zone, U1)] = geom.gCon[0][1]*geom.alpha;
  }

}
