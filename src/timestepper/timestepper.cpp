#include "timestepper.hpp"

timeStepper::timeStepper(const int N1, const int N2, const int N3,
                         const int numGhost, const int dim, 
                         const int numVars, 
                         const double time,
                         const double dt
                        )
{
  this->time = time;
  this->dt = dt;
  this->N1 = N1;
  this->N2 = N2;
  this->N3 = N3;
  this->numVars = numVars;

  XCoords = new grid(N1, N2, N3,
                     numGhost, dim, 3,
                     DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                     DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                     DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                    );

  prim         = new grid(N1, N2, N3,
                          numGhost, dim, numVars,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                         );

  primHalfStep = new grid(N1, N2, N3,
                          numGhost, dim, numVars,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                         );

  primOld      = new grid(N1, N2, N3,
                          numGhost, dim, numVars,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                         );

  cons         = new grid(N1, N2, N3,
                          numGhost, dim, numVars,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                         );

  consOld      = new grid(N1, N2, N3,
                          numGhost, dim, numVars,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                          DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                         );

  sourcesExplicit    = new grid(N1, N2, N3,
                                numGhost, dim, numVars,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                               );

  sourcesImplicit    = new grid(N1, N2, N3,
                                numGhost, dim, numVars,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                               );

  sourcesImplicitOld = new grid(N1, N2, N3,
                                numGhost, dim, numVars,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                               );

  sourcesTimeDer     = new grid(N1, N2, N3,
                                numGhost, dim, numVars,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                               );

  primLeft  = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  primRight = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  fluxesX1  = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  fluxesX2  = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  fluxesX3  = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  divFluxes = new grid(N1, N2, N3,
                       numGhost, dim, numVars,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                       DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                      );

  setXCoords(locations::LEFT,   *XCoords);
  geomLeft    = new geometry(*XCoords);

  setXCoords(locations::RIGHT,  *XCoords);
  geomRight   = new geometry(*XCoords);

  setXCoords(locations::BOTTOM, *XCoords);
  geomBottom  = new geometry(*XCoords);

  setXCoords(locations::TOP,    *XCoords);
  geomTop     = new geometry(*XCoords);

  setXCoords(locations::CENTER, *XCoords);
  geomCenter  = new geometry(*XCoords);
  geomCenter->computeConnectionCoeffs();
  /* XCoords set to locations::CENTER */


  int numReads, numWrites;
  elem          = new fluidElement(prim->vars, *geomCenter,
                                   numReads, numWrites
                                  ); /* n+1   */
  elemOld       = new fluidElement(primOld->vars, *geomCenter,
                                   numReads, numWrites
                                  ); /* n     */
  elemHalfStep  = new fluidElement(primHalfStep->vars, *geomCenter,
                                   numReads, numWrites
                                  ); /* n+1/2 */

  riemann = new riemannSolver(*prim, *geomCenter);

  /* Data structures needed for the nonlinear solver */
  residual        = new grid(N1, N2, N3,
                             numGhost, dim, numVars,
                             DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                             DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                             DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                            );

  residualPlusEps  = new grid(N1, N2, N3,
                              numGhost, dim, numVars,
                              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                             );

  primGuessPlusEps = new grid(N1, N2, N3,
                              numGhost, dim, numVars,
                              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                             );

  primGuessLineSearchTrial = new grid(N1, N2, N3,
                                      numGhost, dim, numVars,
                                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                                     );

  /* The grid data structure arranges data in Struct of Arrays format. Need to
   * rearrange to Array of Structs format in order to solve the linear system Ax
   * = b */
  array zero = 0.*residual->vars[0];

  /* Arrays are copy-on-write. Hence, can use residual->varsSoA to initialize */
  residualSoA  = residual->varsSoA;

  /* Jacobian \partial residual/ \prim in Struct of Arrays format */
  jacobianSoA  = af::constant(1., residual->vars[0].dims(0),
                                  residual->vars[0].dims(1),
                                  residual->vars[0].dims(2),
                                  numVars * numVars,
                                  f64
                             );

  /* Correction dP_k in P_{k+1} = P_k + lambda*dP_k in Array of Structs format */
  deltaPrimAoS = af::constant(0., 
                              numVars,
                              residual->vars[0].dims(0),
                              residual->vars[0].dims(1),
                              residual->vars[0].dims(2),
                              f64
                             );

  /* Steplength lambda in P_{k+1} = P_k + lambda*dP_k */ 
  stepLength   = zero;
                     
  int N1Total = residual->N1Total;
  int N2Total = residual->N1Total;
  int N3Total = residual->N1Total;

  AHostPtr = new double [numVars*numVars*N1Total*N2Total*N3Total];
  bHostPtr = new double [numVars*N1Total*N2Total*N3Total];
  xHostPtr = new double [numVars*N1Total*N2Total*N3Total];

  initialConditions(XCoords->vars, primOld->vars);
}

timeStepper::~timeStepper()
{
  delete prim, primHalfStep, primOld;
  delete cons, consOld;
  delete sourcesExplicit, sourcesImplicit, sourcesImplicitOld, sourcesTimeDer;
  delete primLeft, primRight;
  delete fluxesX1, fluxesX2, fluxesX3;
  delete elem, elemOld, elemHalfStep;
  delete riemann;
  delete geomLeft, geomRight, geomBottom, geomTop, geomCenter;

  delete primGuessLineSearchTrial;
  delete primGuessPlusEps;
  delete residual;
  delete residualPlusEps;

  delete[] AHostPtr;
  delete[] bHostPtr;
  delete[] xHostPtr;
}
