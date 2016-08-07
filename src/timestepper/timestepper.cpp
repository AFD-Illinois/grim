#include "timestepper.hpp"

timeStepper::timeStepper(const int N1, 
                         const int N2,
                         const int N3,
                         const int dim,
                         const int numVars, 
                         const int numGhost,
                         const double time,
                         const double dt,
                         const int boundaryLeft,  const int boundaryRight,
                         const int boundaryTop,   const int boundaryBottom,
                         const int boundaryFront, const int boundaryBack,
                         const int metric,
                         const double blackHoleSpin,
                         const double hSlope,
                         const double X1Start, const double X1End,
                         const double X2Start, const double X2End,
                         const double X3Start, const double X3End
                        )
{
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);

  PetscPrintf(PETSC_COMM_WORLD, "   _____ _____  _____ __  __ \n");
  PetscPrintf(PETSC_COMM_WORLD, "  / ____|  __ \\|_   _|  \\/  |\n");
  PetscPrintf(PETSC_COMM_WORLD,  " | |  __| |__) | | | | \\  / |\n");
  PetscPrintf(PETSC_COMM_WORLD,  " | | |_ |  _  /  | | | |\\/| |\n");
  PetscPrintf(PETSC_COMM_WORLD,  " | |__| | | \\ \\ _| |_| |  | |\n");
  PetscPrintf(PETSC_COMM_WORLD,  "  \\_____|_|  \\_\\_____|_|  |_|\n");
  PetscPrintf(PETSC_COMM_WORLD,  "\n");
  
  this->time = time;
  this->dt = dt;
  this->numGhost = numGhost;
  this->dim = dim;
  this->numVars = numVars;

  this->boundaryLeft   = boundaryLeft;
  this->boundaryRight  = boundaryRight;
  this->boundaryTop    = boundaryTop;
  this->boundaryBottom = boundaryBottom;
  this->boundaryFront  = boundaryFront;
  this->boundaryBack   = boundaryBack;

  int periodicBoundariesX1;
  int periodicBoundariesX2;
  int periodicBoundariesX3;

  if (   boundaryLeft  == boundaries::PERIODIC
      || boundaryRight == boundaries::PERIODIC
     )
  {
    periodicBoundariesX1 = 1;
  }

  if (   boundaryTop    == boundaries::PERIODIC
      || boundaryBottom == boundaries::PERIODIC
     )
  {
    periodicBoundariesX2 = 1;
  }

  if (   boundaryFront == boundaries::PERIODIC
      || boundaryBack  == boundaries::PERIODIC
     )
  {
    periodicBoundariesX3 = 1;
  }

  prim         = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  /* Set N1, N2, N3 here after creating a grid(), so that they are correctly set
   * according to dim */
  this->N1 = prim->N1;
  this->N2 = prim->N2;
  this->N3 = prim->N3;
  std::string deviceInfo = af::infoString();
  
  double availableBandwidth = bandwidthTest(10000);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, 
  "#### Rank %d of %d: System info ####\n %s \n  Local size       : %i x %i x %i\n  Memory Bandwidth : %g GB/sec\n\n",
                          world_rank, world_size, deviceInfo.c_str(),
                          prim->N1Total, prim->N2Total, prim->N3Total,
                          availableBandwidth
                         );  
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);


  primHalfStep = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  primOld      = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  cons         = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  consOld      = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  sourcesExplicit    = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  sourcesImplicit    = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  sourcesImplicitOld = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  sourcesTimeDer     = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  primLeft  = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  primRight = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  fluxesX1  = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  fluxesX2  = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  fluxesX3  = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  emfX1  = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  emfX2  = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  emfX3  = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  divFluxes = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  divB = new grid(N1, N2, N3,
                  dim, 1, numGhost,
                  periodicBoundariesX1,
                  periodicBoundariesX2,
                  periodicBoundariesX3
                 );

  XCoords = new coordinatesGrid(N1, N2, N3,
                                dim, numGhost,
                                X1Start, X1End,
                                X2Start, X2End,
                                X3Start, X3End
                               );

  PetscPrintf(PETSC_COMM_WORLD, "  Generating metric at LEFT face...");
  XCoords->setXCoords(locations::LEFT);
  geomLeft    = new geometry(metric,
                             blackHoleSpin,
                             hSlope, 
                             *XCoords
                            );
  PetscPrintf(PETSC_COMM_WORLD, "done\n");

  PetscPrintf(PETSC_COMM_WORLD, "  Generating metric at RIGHT face...");
  XCoords->setXCoords(locations::RIGHT);
  geomRight   = new geometry(metric,
                             blackHoleSpin,
                             hSlope,
                             *XCoords
                            );
  PetscPrintf(PETSC_COMM_WORLD, "done\n");

  PetscPrintf(PETSC_COMM_WORLD, "  Generating metric at BOTTOM face...");
  XCoords->setXCoords(locations::BOTTOM);
  geomBottom  = new geometry(metric,
                             blackHoleSpin,
                             hSlope,
                             *XCoords
                            );
  PetscPrintf(PETSC_COMM_WORLD, "done\n");

  PetscPrintf(PETSC_COMM_WORLD, "  Generating metric at TOP face...");
  XCoords->setXCoords(locations::TOP);
  geomTop     = new geometry(metric,
                             blackHoleSpin,
                             hSlope, 
                             *XCoords
                            );
  PetscPrintf(PETSC_COMM_WORLD, "done\n");

  PetscPrintf(PETSC_COMM_WORLD, "  Generating metric at zone CENTER...");
  XCoords->setXCoords(locations::CENTER);
  geomCenter  = new geometry(metric,
                             blackHoleSpin,
                             hSlope,
                             *XCoords
                            );
  PetscPrintf(PETSC_COMM_WORLD, "done\n");

  PetscPrintf(PETSC_COMM_WORLD, "  Computing connections at zone CENTER...");
  geomCenter->computeConnectionCoeffs();
  PetscPrintf(PETSC_COMM_WORLD, "done\n\n");
  /* XCoords set to locations::CENTER */

  int numReads, numWrites;
  elem          = new fluidElement(*prim, *geomCenter,
                                   numReads, numWrites
                                  ); /* n+1   */
  elemOld       = new fluidElement(*primOld, *geomCenter,
                                   numReads, numWrites
                                  ); /* n     */
  elemHalfStep  = new fluidElement(*primHalfStep, *geomCenter,
                                   numReads, numWrites
                                  ); /* n+1/2 */

  riemann = new riemannSolver(*prim, *geomCenter);

  int numFluidVars = vars::numFluidVars;
  /* Data structures needed for the nonlinear solver */
  residual        = new grid(N1, N2, N3,
                             dim, numFluidVars, numGhost,
                             periodicBoundariesX1,
                             periodicBoundariesX2,
                             periodicBoundariesX3
                            );

  residualPlusEps  = new grid(N1, N2, N3,
                              dim, numFluidVars, numGhost,
                              periodicBoundariesX1,
                              periodicBoundariesX2,
                              periodicBoundariesX3
                             );

  primGuessPlusEps = new grid(N1, N2, N3,
                              dim, numVars, numGhost,
                              periodicBoundariesX1,
                              periodicBoundariesX2,
                              periodicBoundariesX3
                             );

  primGuessLineSearchTrial = new grid(N1, N2, N3,
                                      dim, numVars, numGhost,
                                      periodicBoundariesX1,
                                      periodicBoundariesX2,
                                      periodicBoundariesX3
                                     );

  primIC = new grid(N1, N2, N3,
		    dim, numVars, numGhost,
		    periodicBoundariesX1,
		    periodicBoundariesX2,
		    periodicBoundariesX3
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
                                  numFluidVars * numFluidVars,
                                  f64
                             );

  /* Correction dP_k in P_{k+1} = P_k + lambda*dP_k in Array of Structs format */
  deltaPrimAoS = af::constant(0., 
                              numFluidVars,
                              residual->vars[0].dims(0),
                              residual->vars[0].dims(1),
                              residual->vars[0].dims(2),
                              f64
                             );

  /* Steplength lambda in P_{k+1} = P_k + lambda*dP_k */ 
  stepLength  = zero;
                     
  int N1Total = residual->N1Total;
  int N2Total = residual->N2Total;
  int N3Total = residual->N3Total;

  AHostPtr = new double [numFluidVars*numFluidVars*N1Total*N2Total*N3Total];
  bHostPtr = new double [numFluidVars*N1Total*N2Total*N3Total];

  /* Mask for ghost zone residuals */
  residualMask = af::constant(0.,
                              residual->vars[0].dims(0),
                              residual->vars[0].dims(1),
                              residual->vars[0].dims(2),
                              f64
                             );
  /* Get the domain of the bulk */
  domainX1 = *residual->domainX1;
  domainX2 = *residual->domainX2;
  domainX3 = *residual->domainX3;

  residualMask(domainX1, domainX2, domainX3) = 1.;

  initialConditions(numReads, numWrites);
  if (params::restart)
  {
    struct stat fileInfo;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if (rank==0)
    {
      if (stat(params::restartFile.c_str(), &fileInfo) == 0)
      {
        /* File exists */
        PetscPrintf(PETSC_COMM_WORLD, "\nFound restart file: %s\n\n", 
                    params::restartFile.c_str()
                   ); 
      }
      else
      {
        /* File does not exist */
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        PetscPrintf(PETSC_COMM_WORLD, "Restart file %s does not exist\n",
                    params::restartFile.c_str()
                   );
        MPI_Abort(PETSC_COMM_WORLD, 1);
      }
    }

    primOld->load("primitives", params::restartFile);
    // Need to call the diagnostics to reset the time step !!!
    fullStepDiagnostics(numReads, numWrites);
  }
}

timeStepper::~timeStepper()
{
  delete prim, primHalfStep, primOld;
  delete cons, consOld;
  delete sourcesExplicit, sourcesImplicit, sourcesImplicitOld, sourcesTimeDer;
  delete primLeft, primRight;
  delete fluxesX1, fluxesX2, fluxesX3;
  delete divFluxes;
  delete divB;
  delete emfX1, emfX2, emfX3;
  delete elem, elemOld, elemHalfStep;
  delete riemann;
  delete geomLeft, geomRight, geomBottom, geomTop, geomCenter;

  delete primGuessLineSearchTrial;
  delete primGuessPlusEps;
  delete residual;
  delete residualPlusEps;

  delete[] AHostPtr;
  delete[] bHostPtr;
}

/* Returns memory bandwidth in GB/sec */
double timeStepper::memoryBandwidth(const double numReads,
                                    const double numWrites,
                                    const double numEvals,
                                    const double timeElapsed
                                    )
{
  return   (double)(prim->N1Total)
         * (double)(prim->N2Total)
         * (double)(prim->N3Total)
         * 8. * (numReads + numWrites) * numEvals /timeElapsed/1e9;
}

double timeStepper::bandwidthTest(const int numEvals)
{
  grid prim(N1, N2, N3,
            dim, 3, numGhost,
            0, 0, 0
           );

  prim.vars[0] = af::randu(prim.vars[0].dims(0),
                           prim.vars[0].dims(1),
                           prim.vars[0].dims(2),
                           f64
                          );
  prim.vars[1] = af::randu(prim.vars[1].dims(0),
                           prim.vars[1].dims(1),
                           prim.vars[1].dims(2),
                           f64
                          );
  prim.vars[2] = af::randu(prim.vars[2].dims(0),
                           prim.vars[2].dims(1),
                           prim.vars[2].dims(2),
                           f64
                          );
  prim.vars[2] = prim.vars[1] + prim.vars[0];
  prim.vars[2].eval();
  af::sync();

  af::timer::start();
  for (int n=0; n < numEvals; n++)
  {
    prim.vars[2] = prim.vars[1] + prim.vars[0];
    prim.vars[2].eval();
  }
  af::sync();

  int numReads  = 2;
  int numWrites = 1;
  double timeElapsed = af::timer::stop();

  return memoryBandwidth(numReads, numWrites,
                         numEvals, timeElapsed
                        );
}

