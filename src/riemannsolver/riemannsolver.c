#include "riemannsolver.h"

void riemannSolver(const REAL fluxLeft[ARRAY_ARGS DOF],
                   const REAL fluxRight[ARRAY_ARGS DOF],
                   const REAL conservedVarsLeft[ARRAY_ARGS DOF],
                   const REAL conservedVarsRight[ARRAY_ARGS DOF],
                   const REAL primVarsLeft[ARRAY_ARGS DOF],
                   const REAL primVarsRight[ARRAY_ARGS DOF],
                   const struct geometry geom[ARRAY_ARGS DOF],
                   const int dir, REAL fluxes[ARRAY_ARGS DOF])
{
  REAL cMinLeft, cMaxLeft;
  struct fluidElement elem;
  setFluidElement(primVarsLeft, geom, &elem);
  waveSpeeds(&elem, geom, dir, &cMinLeft, &cMaxLeft);

  REAL cMinRight, cMaxRight;
  setFluidElement(primVarsRight, geom, &elem);
  waveSpeeds(&elem, geom, dir, &cMinRight, &cMaxRight);
  
	REAL cMax = fabs(fmax(fmax(0., cMaxLeft), cMaxRight));
	REAL cMin = fabs(fmax(fmax(0., -cMinLeft), -cMinRight));
	REAL cLaxFriedrichs = fmax(cMax, cMin);
    
  for (int var=0; var<DOF; var++) 
  {
    fluxes[var] = 0.5*(fluxLeft[var] + fluxRight[var]
                       - cLaxFriedrichs*(  conservedVarsRight[var]
                                         - conservedVarsLeft[var]
                                        )
                    );
                         
  }
}

void waveSpeeds(const struct fluidElement elem[ARRAY_ARGS 1],
                const struct geometry geom[ARRAY_ARGS 1],
                const int dir,
                REAL cMin[ARRAY_ARGS 1], REAL cMax[ARRAY_ARGS 1])
{
  REAL bCov[NDIM];
  conToCov(elem->bCon, geom, bCov);

  REAL bSqr = covDotCon(bCov, elem->bCon);

  REAL cAlvenSqr = bSqr/(bSqr + elem->primVars[RHO] 
                              + ADIABATIC_INDEX*elem->primVars[UU]);
  REAL csSqr = (ADIABATIC_INDEX)*(ADIABATIC_INDEX-1)*elem->primVars[UU]
              /(elem->primVars[RHO] + ADIABATIC_INDEX*elem->primVars[UU]);

  REAL cmSqr = csSqr + cAlvenSqr - csSqr*cAlvenSqr;
  
  REAL ACov[NDIM], ACon[NDIM];
  REAL BCov[NDIM], BCon[NDIM];
  for (int mu=0; mu<NDIM; mu++)
  {
    ACov[mu] = 0.;
    BCov[mu] = 0.;
  }
  ACov[dir] = 1.;
  BCov[0] = 1.;
  covToCon(ACov, geom, ACon);
  covToCon(BCov, geom, BCon);

  REAL ASqr = 0., BSqr=0., ADotU=0., BDotU=0., ADotB=0.; 
  for (int mu=0; mu<NDIM; mu++)
  {
    ASqr += ACov[mu]*ACon[mu];
    BSqr += BCov[mu]*BCon[mu];
    ADotU += ACov[mu]*elem->uCon[mu];
    BDotU += BCov[mu]*elem->uCon[mu];
    ADotB += ACov[mu]*BCon[mu];
  }

  REAL A = (BDotU*BDotU) - (BSqr + BDotU*BDotU)*cmSqr;
  REAL B = 2.*(ADotU*BDotU - (ADotB + ADotU*BDotU)*cmSqr);
  REAL C = ADotU*ADotU - (ASqr + ADotU*ADotU)*cmSqr;

  REAL discr = sqrt(B*B - 4.*A*C);

  REAL cPlus = -(-B + discr)/(2.*A);
  REAL cMinus = -(-B - discr)/(2.*A);

  if (cPlus>cMinus) {
    cMax[0] = cPlus;
    cMin[0] = cMinus;
  } else {
    cMax[0] = cMinus;
    cMin[0] = cPlus;
  }

}
