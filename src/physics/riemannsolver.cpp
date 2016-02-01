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

  fluxLeft = new grid(N1, N2, N3,
                      numGhost, dim, numVars,
                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                     );

  fluxRight = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  consLeft = new grid(N1, N2, N3,
                      numGhost, dim, numVars,
                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                     );

  consRight = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  int numReads, numWrites;
  elemFace  = new fluidElement(prim.vars, geom,
                               numReads, numWrites
                              );

//  MinSpeedLeft = af::constant(1,
//			  primLeft->vars[0].dims(directions::X1),
//			  primLeft->vars[0].dims(directions::X2),
//			  primLeft->vars[0].dims(directions::X3),
//			  f64
//			  );
//  MaxSpeedLeft = af::constant(1,
//			  primLeft->vars[0].dims(directions::X1),
//			  primLeft->vars[0].dims(directions::X2),
//			  primLeft->vars[0].dims(directions::X3),
//			  f64
//			  );
//  MinSpeedRight = af::constant(1,
//			  primLeft->vars[0].dims(directions::X1),
//			  primLeft->vars[0].dims(directions::X2),
//			  primLeft->vars[0].dims(directions::X3),
//			  f64
//			  );
//  MaxSpeedRight = af::constant(1,
//			  primLeft->vars[0].dims(directions::X1),
//			  primLeft->vars[0].dims(directions::X2),
//			  primLeft->vars[0].dims(directions::X3),
//			  f64
//			  );
}

riemannSolver::~riemannSolver()
{
  delete fluxLeft, fluxRight;
  delete consLeft, consRight;
  delete elemFace;
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

//  /* Reconstruction gives, at a point of index i:
//   * primLeft : right-biased stencil reconstructs on face i-/1.2
//   * primRight: left-biased stencil reconstructs on face i+1/2 */
//  reconstruction::reconstruct(prim, dir, *primLeft, *primRight);

  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;

  /* Compute fluxes and cons at i+1/2 - eps : left flux on right face */
  elemFace->set(primRight.vars, geomRight,
                numReadsElemSet, numWritesElemSet
               );
  elemFace->computeFluxes(geomRight, fluxDirection, fluxLeft->vars,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );
  elemFace->computeFluxes(geomRight, 0,             consLeft->vars,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );

  /* Compute fluxes and cons at i-1/2 + eps : right flux on left face */
  elemFace->set(primLeft.vars, geomLeft,
                numReadsElemSet, numWritesElemSet
               );
  elemFace->computeFluxes(geomLeft, fluxDirection, fluxRight->vars,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );
  elemFace->computeFluxes(geomLeft, 0,             consRight->vars,
                          numReadsComputeFluxes, numWritesComputeFluxes
                         );

//  elemLeft->computeMinMaxCharSpeeds(geom,dir,MinSpeedLeft,MaxSpeedLeft);
//  elemRight->computeMinMaxCharSpeeds(geom,dir,MinSpeedRight,MaxSpeedRight);
  
  numReads = 2*(numReadsElemSet + 2*numReadsComputeFluxes);
  numWrites = 2*(numWritesElemSet + 2*numWritesComputeFluxes);

  for (int var=0; var<primLeft.numVars; var++)
  {
    //The fluxes are requested on the left-face i-1/2.
    //Hence, we need to shift elemLeft,fluxLeft,consLeft by a single point to the right
    // (elemLeft[i] refers to values at i+1/2, we want values at i-1/2)
//    MinSpeedLeft = af::shift(MinSpeedLeft,shiftX1, shiftX2, shiftX3);
//    MaxSpeedLeft = af::shift(MaxSpeedLeft,shiftX1, shiftX2, shiftX3);
//    MinSpeedLeft = min(MinSpeedLeft,MinSpeedRight);
//    MaxSpeedLeft = max(MaxSpeedLeft,MaxSpeedRight);

    //HLL formula
//    flux.vars[var] = 
//      (  maxSpeedLeft * af::shift(fluxLeft->vars[var], shiftX1, shiftX2, shiftX3) 
//			 - minSpeedLeft * fluxRight->vars[var]
//			 + minSpeedLeft * maxSpeedLeft 
//			 * (  consRight->vars[var]
//			     - af::shift(consLeft->vars[var], shiftX1, shiftX2, shiftX3)
//			   )
//			)/(maxSpeedLeft - minSpeedLeft);

    flux.vars[var] = 
      af::shift(fluxLeft->vars[var], shiftX1, shiftX2, shiftX3)
    - fluxRight->vars[var]
    + consRight->vars[var] 
    - af::shift(consLeft->vars[var], shiftX1, shiftX2, shiftX3);

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
