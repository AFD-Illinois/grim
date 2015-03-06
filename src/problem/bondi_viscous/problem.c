#include "../problem.h"
#include "bondi_viscous.h"

REAL FuncT(const REAL T, const REAL R, const REAL C1, const REAL C2)
{
    REAL nPoly = 1./(ADIABATIC_INDEX-1.);
    return pow(1.+(1.+nPoly)*T,2.)*(1.-2./R+pow(C1/R/R/pow(T,nPoly),2.))-C2;
}

REAL SolveT(const REAL R, const REAL C1, const REAL C2)
{
  REAL nPoly = 1./(ADIABATIC_INDEX-1.);
  REAL rtol = 1.e-12;
  REAL ftol = 1.e-14;
  REAL Tmin = 0.6*(sqrt(C2)-1.)/(nPoly+1.);
  REAL Tmax = pow(C1*sqrt(2./R/R/R),1./nPoly);
  REAL f0,f1,fh;
  REAL T0,T1,Th;
  T0=0.6*Tmin;
  f0=FuncT(T0,R,C1,C2);
  T1=Tmax;
  f1=FuncT(T1,R,C1,C2);
  if(f0*f1>0.)
    {
      printf("Failed solving for T at R = %f; C1 = %f; C2 = %f \n",R,C1,C2);
      PetscPrintf(PETSC_COMM_WORLD, "Failed determination of T \n");
      exit(1);
    }
  Th = (f1*T0-f0*T1)/(f1-f0);
  fh = FuncT(Th,R,C1,C2);
  REAL EpsT = rtol*(Tmin+Tmax);
  while(fabs(Th-T0)>EpsT && fabs(Th-T1)>EpsT && fabs(fh)>ftol)
    {
      if(fh*f0<0.)
	{
	  T0=Th;
	  f0=fh;
	}
      else
	{
	  T1=Th;
	  f1=fh;
	}
      Th = (f1*T0-f0*T1)/(f1-f0);
      fh = FuncT(Th,R,C1,C2);
    }
  return Th;
}


void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  PetscPrintf(PETSC_COMM_WORLD, "Initializing hydro variables...");

  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  //Constanst of the Bondi solution from location of sonic point,
  //following Hawley, Smarr & Wilson 1984
  REAL nPoly = 1./(ADIABATIC_INDEX-1.); 
  REAL Rc = SONICRADIUS;
  REAL uc = sqrt(0.5/Rc);
  REAL vc2 = nPoly*uc*uc/(1.-3.*uc*uc);
  REAL Tc = vc2/(1.-vc2)/(nPoly+1.);
  REAL C1 = pow(Tc,nPoly)*uc*Rc*Rc;
  REAL C2 = pow(1.+(nPoly+1.)*Tc,2.)*(1.-2./Rc+uc*uc);
  REAL Kp = pow(MDOT/4./M_PI/C1,-1./nPoly);

  REAL randNum;
  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
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

      REAL r = xCoords[1];
      REAL T = 1.;
      REAL Rho0 = RHO_FLOOR_MIN;
      REAL ur = 0.;
      REAL ut = 1.;
      if(r>2.+1.e-8)
	{
	  T = SolveT(r,C1,C2);
	  Rho0 = pow(T/Kp,nPoly);
	  ur = -C1/r/r/pow(T,nPoly);
	  ut = (2./r*ur+sqrt(1.-2./r+ur*ur))/(1.-2./r);
	}

      PetscRandomGetValue(randNumGen, &randNum);

      primTile[INDEX_TILE(&zone, RHO)] = Rho0;
      primTile[INDEX_TILE(&zone, UU)] = nPoly*Rho0*T;
      primTile[INDEX_TILE(&zone, U1)] = ur + ut*2./r+.001*(randNum-0.5);
      primTile[INDEX_TILE(&zone, U2)] = 0.;
      primTile[INDEX_TILE(&zone, U3)] = 0.;

      primTile[INDEX_TILE(&zone, B1)] = BMAG/r/r/sqrt(1.+2./r);;
      primTile[INDEX_TILE(&zone, B2)] = 0.;
      primTile[INDEX_TILE(&zone, B3)] = 0.;

      //Dirichlet Boundary data
      if (zone.i < 0)
      {
        for (int var=0; var<DOF; var++)
        {
          ts->problemSpecificData->primVarsLeftEdge[zone.i+NG][var]
            = primTile[INDEX_TILE(&zone, var)];
        }
      }
      if (zone.i > N1-1)
      {
        for (int var=0; var<DOF; var++)
        {
          ts->problemSpecificData->primVarsRightEdge[zone.i-N1][var]
            = primTile[INDEX_TILE(&zone, var)];
        }
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
      
      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primOldGlobal, &zone, var) =
          primTile[INDEX_TILE(&zone, var)];
      }
    }

  }
  PetscRandomDestroy(&randNumGen);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);

  PetscPrintf(PETSC_COMM_WORLD, "done\n");
  /* Done with setting the initial conditions */

}

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
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

    if (primTile[INDEX_TILE(&zone, RHO)] < RHO_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR_MIN;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < UU_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR_MIN;
    }
  }
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

    REAL rho = primTile[INDEX_TILE(&zone, RHO)];
    REAL u = primTile[INDEX_TILE(&zone, UU)];

    if (rho < rhoFloor)
    {
      primTile[INDEX_TILE(&zone, RHO)] = rhoFloor;
    }

    if (u < uFloor)
    {
      primTile[INDEX_TILE(&zone, UU)] = uFloor;
    }

    struct geometry geom;
    setGeometry(XCoords, &geom);

    struct fluidElement elem;
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);

    if (elem.gamma > GAMMA_MAX)
    {
      REAL factor = sqrt( (GAMMA_MAX*GAMMA_MAX-1.)
                         /(elem.gamma*elem.gamma-1.)
                        );

      primTile[INDEX_TILE(&zone, U1)] *= factor;
      primTile[INDEX_TILE(&zone, U2)] *= factor;
      primTile[INDEX_TILE(&zone, U3)] *= factor;
    }
  }

}

void applyAdditionalProblemSpecificBCs
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL primTile[ARRAY_ARGS TILE_SIZE]
)
{
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    #if (PHYSICAL_BOUNDARY_LEFT_EDGE==DIRICHLET)
      if (zone.i < 0)
      {
        for (int var=0; var<DOF; var++)
        {
          primTile[INDEX_TILE(&zone, var)] =
            problemSpecificData->primVarsLeftEdge[zone.i+NG][var];
        }

      }
    #endif

    #if (PHYSICAL_BOUNDARY_RIGHT_EDGE==DIRICHLET)
      if (zone.i > N1-1)
      {
        for (int var=0; var<DOF; var++)
        {
          primTile[INDEX_TILE(&zone, var)] = 
            problemSpecificData->primVarsRightEdge[zone.i-N1][var];
        }
      }
    #endif

  }

}

void applyProblemSpecificFluxFilter
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
)
{
}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
}

void inflowCheck(const struct gridZone zone[ARRAY_ARGS 1],
                 REAL primTile[ARRAY_ARGS TILE_SIZE])
{
}

#if (CONDUCTION)
void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
                             struct fluidElement elem[ARRAY_ARGS 1])
{
  SETERRQ(PETSC_COMM_WORLD, 1,
          "Conduction parameters not set in shock_tests/problem.c\n");
}
#endif
