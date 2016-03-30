#include "physics.hpp"

riemannSolver::riemannSolver(const grid &prim,
                             const geometry &geom
                            )
{
  int N1 = prim.N1;
  int N2 = prim.N2;
  int N3 = prim.N3;

  int numGhost = prim.numGhost;
  int dim      = prim.dim;
  int numVars  = prim.numVars;

  fluxLeft  = new grid(N1, N2, N3, 
                       dim, numVars, numGhost,
                       false, false, false
                      );

  fluxRight = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       false, false, false
                      );

  consLeft  = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       false, false, false
                      );

  consRight = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       false, false, false
                      );

  int numReads, numWrites;
  elemFace  = new fluidElement(prim, geom,
                               numReads, numWrites
                              );

  /* Allocate space for the wavespeeds using elemFace->one */
  minSpeedLeft  = elemFace->one;
  minSpeedRight = elemFace->one;
  maxSpeedLeft  = elemFace->one;
  maxSpeedRight = elemFace->one;
}

riemannSolver::~riemannSolver()
{
  delete fluxLeft, fluxRight;
  delete consLeft, consRight;
  delete elemFace;
}

void fluidElement::computeMinMaxCharSpeeds(const geometry &geom,
			                                     const int dir,
                                			     array &minSpeed,
                                			     array &maxSpeed,
                                           int &numReads,
                                           int &numWrites
                                          )
{
  array zero = 0.*one;
  numReads = 0;
  numWrites = 0;
  
  /* Alven speed */
  array cAlvenSqr = bSqr/(rho+params::adiabaticIndex*u+bSqr);
  numReads += 3;
  /* Reads:
   * -----
   *  bSqr, rho, u : 3
   *
   * Writes: 0
   * ------ */

  /* Hydro sound speed */
  array csSqr = soundSpeed*soundSpeed;
  numReads += 1;
  /* Reads:
   * -----
   *  soundSpeed : 1
   *
   * Writes: 0
   * ------ */

  /* Approximate contribution from dP */
  array cVisSqr = zero;
  if (params::viscosity)
  {
    cVisSqr = 4./3./(rho+params::adiabaticIndex*u)*rho*nu_emhd/tau;
    numReads += 2;
    /* Reads:
     * -----
     *  nu_emhd, tau : 2
     *
     * Writes: 0
     * ------ */
  }

  /* Approximate contribution from Q */
  array cConSqr = zero;
  if (params::conduction)
  {
    cConSqr = (params::adiabaticIndex-1.)*chi_emhd/tau;
    numReads += 1;
    /* Reads:
     * -----
     *  nu_emhd : 1
     *  (Not counting tau since it is already read if params::viscosity==1)
     *
     * Writes: 0
     * ------ */
  }
  cConSqr = 0.5*(csSqr+cConSqr+af::sqrt(csSqr*csSqr+cConSqr*cConSqr));
  csSqr = cConSqr + cVisSqr;

  /* Combine speeds using relativistic velocity additions */
  csSqr = csSqr + cAlvenSqr - csSqr*cAlvenSqr; 
  
  array condition = csSqr > one;
  csSqr = csSqr*(one-condition)+condition;
  
  int sdir = 0;
  switch(dir)
  {
    case directions::X1:
      sdir=1; break;
    case directions::X2:
      sdir=2; break;
    case directions::X3:
      sdir=3; break;
  }

  array ACov[NDIM], ACon[NDIM];
  array BCov[NDIM], BCon[NDIM];
  for (int mu=0; mu<NDIM; mu++)
  {
    ACov[mu] = zero;
    BCov[mu] = zero;
  }
  ACov[sdir] = one;
  BCov[0]    = one;
  for (int mu=0; mu<NDIM; mu++)
  {
    ACon[mu] = zero;
    BCon[mu] = zero;

    for(int nu=0;nu<NDIM; nu++)
	  {
	    ACon[mu] += geom.gCon[mu][nu]*ACov[nu];
	    BCon[mu] += geom.gCon[mu][nu]*BCov[nu];
	  }
  }
  numReads += 16;
  /* Reads:
   * -----
   * geom.gCon[mu][nu] : 16
   *
   * Writes: 0
   * ------ */

  array ASqr = zero;
  array BSqr = zero;
  array ADotU = zero;
  array BDotU = zero;
  array ADotB = zero;
  for (int mu=0; mu<NDIM; mu++)
  {
    ASqr += ACov[mu]*ACon[mu];
    BSqr += BCov[mu]*BCon[mu];
    ADotU += ACov[mu]*uCon[mu];
    BDotU += BCov[mu]*uCon[mu];
    ADotB += ACov[mu]*BCon[mu];
  }
  
  ASqr.eval();
  BSqr.eval();
  ADotU.eval();
  BDotU.eval();
  ADotB.eval();

  array A = (BDotU*BDotU)   - (BSqr + BDotU*BDotU)*csSqr;
  array B = 2.*(ADotU*BDotU - (ADotB + ADotU*BDotU)*csSqr);
  array C = ADotU*ADotU     - (ASqr + ADotU*ADotU)*csSqr;
  array discr = af::sqrt(B*B - 4.*A*C);

  minSpeed = -(-B + discr)/2./A;
  maxSpeed = -(-B - discr)/2./A;
  
  condition = (minSpeed > -1.e-15);
  minSpeed  = minSpeed*(1-condition)-1.e-15*condition;
  minSpeed.eval();

  condition = (maxSpeed < 1.e-15);
  maxSpeed  = maxSpeed*(1-condition)+1.e-15*condition;
  maxSpeed.eval();
  numWrites += 2;
  /* Reads: 0
   * -----
   *
   * Writes:
   * minSpeed, maxSpeed : 2
   * ------ */
}

void riemannSolver::solve(const grid &primLeft,
                          const grid &primRight,
                          const geometry &geomLeft,
                          const geometry &geomRight,
                          const int dir,
                          grid &flux,
                          int &numReads,
                          int &numWrites
                         )
{
  int shiftX1, shiftX2, shiftX3;
  int fluxDirection;
  switch (dir)
  {
    case directions::X1:
      fluxDirection = 1;
      shiftX1  = 1;
      shiftX2  = 0;
      shiftX3  = 0;
      break;

    case directions::X2:
      fluxDirection = 2;
      shiftX1  = 0;
      shiftX2  = 1;
      shiftX3  = 0;
      break;

    case directions::X3:
      fluxDirection = 3;
      shiftX1  = 0;
      shiftX2  = 0;
      shiftX3  = 1;
      break;
  }

  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;
  int numReadsCharSpeeds, numWritesCharSpeeds;

  /* Compute fluxes and cons at i+1/2 - eps : left flux on right face */
  elemFace->set(primRight, geomRight,
                numReadsElemSet, numWritesElemSet
               );
  elemFace->computeFluxes(geomRight, fluxDirection, *fluxLeft,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );
  elemFace->computeFluxes(geomRight, 0,             *consLeft,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );
  elemFace->computeMinMaxCharSpeeds(geomRight, dir, 
                                    minSpeedLeft, maxSpeedLeft,
                                    numReadsCharSpeeds, numWritesCharSpeeds
                                   );

  /* Compute fluxes and cons at i-1/2 + eps : right flux on left face */
  elemFace->set(primLeft, geomLeft,
                numReadsElemSet, numWritesElemSet
               );
  elemFace->computeFluxes(geomLeft, fluxDirection, *fluxRight,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );
  elemFace->computeFluxes(geomLeft, 0,             *consRight,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );
  elemFace->computeMinMaxCharSpeeds(geomLeft, dir, 
                                    minSpeedRight, maxSpeedRight,
                                    numReadsCharSpeeds, numWritesCharSpeeds
                                   );

  numReads = 2*(  numReadsElemSet + 2*numReadsComputeFluxes 
                + numReadsCharSpeeds
               );
  numWrites = 2*(numWritesElemSet + 2*numWritesComputeFluxes
                 + numWritesCharSpeeds
                );

  /* The fluxes are requested on the left-face i-1/2.
   * Hence, we need to shift elemLeft,fluxLeft,consLeft by a single point to the right
   * (elemLeft[i] refers to values at i+1/2, we want values at i-1/2) */
  minSpeedLeft = af::shift(minSpeedLeft,shiftX1, shiftX2, shiftX3);
  maxSpeedLeft = af::shift(maxSpeedLeft,shiftX1, shiftX2, shiftX3);
  minSpeedLeft = af::min(minSpeedLeft, minSpeedRight);
  maxSpeedLeft = af::max(maxSpeedLeft, maxSpeedRight);
  numReads += 4;
  /* Reads:
   * -----
   *  minSpeedLeft, minSpeedRight, maxSpeedLeft, maxSpeedRight : 4
   *
   * Writes: 0
   * ------ */

  for (int var=0; var < primLeft.numVars; var++)
  {
    if (params::riemannSolver == riemannSolvers::HLL)
    {
      flux.vars[var] = 
        (   maxSpeedLeft * af::shift(fluxLeft->vars[var], shiftX1, shiftX2, shiftX3) 
			    - minSpeedLeft * fluxRight->vars[var]
			    + minSpeedLeft * maxSpeedLeft 
			    * (  consRight->vars[var]
			       - af::shift(consLeft->vars[var], shiftX1, shiftX2, shiftX3)
			      )
			  )/(maxSpeedLeft - minSpeedLeft);
    }
    else if (params::riemannSolver == riemannSolvers::LOCAL_LAX_FRIEDRICH)
    {
      flux.vars[var] = 
        0.5*(af::shift(fluxLeft->vars[var], shiftX1, shiftX2, shiftX3)
	     + fluxRight->vars[var])
	- af::max(maxSpeedLeft,-minSpeedLeft)*
	(consRight->vars[var] 
	 - af::shift(consLeft->vars[var], shiftX1, shiftX2, shiftX3));
    }

    flux.vars[var].eval();
  }
  /* Reads:
   * -----
   *  fluxLeft[var], fluxRight[var], consLeft[var], consRight[var] : 4*numVars
   *
   * Writes:
   * ------
   * flux[var] : numVars */
  numReads  += 4*primLeft.numVars;
  numWrites +=   primLeft.numVars;
}
