#include "timestepper.hpp"

void timeStepper::timeStep(int &numReads, int &numWrites)
{
  af::timer timeStepTimer = af::timer::start();
  PetscPrintf(PETSC_COMM_WORLD, "  Time = %f, dt = %f\n\n", time, dt);
  int numReadsDt, numWritesDt;
  computeDt(numReadsDt, numWritesDt);

  /* First take a half step */
  PetscPrintf(PETSC_COMM_WORLD, "  ---Half step--- \n");
  af::timer halfStepTimer = af::timer::start();

  currentStep = timeStepperSwitches::HALF_STEP;
  /* Apply boundary conditions on primOld */
  af::timer boundaryTimer = af::timer::start();
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *primOld
                                     );
  setProblemSpecificBCs(numReads,numWrites);
  double boundaryTime = af::timer::stop(boundaryTimer);

  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;
  af::timer elemOldTimer = af::timer::start();
  elemOld->set(*primOld, *geomCenter,
               numReadsElemSet, numWritesElemSet
              );
  double elemOldTime = af::timer::stop(elemOldTimer);

  af::timer consOldTimer = af::timer::start();
  elemOld->computeFluxes(*geomCenter, 0, *consOld, 
                         numReadsComputeFluxes, numWritesComputeFluxes
                        );
  numReads  = numReadsElemSet  + numReadsComputeFluxes;
  numWrites = numWritesElemSet + numWritesComputeFluxes; 
  double consOldTime = af::timer::stop(consOldTimer);

  int numReadsExplicitSouces, numWritesExplicitSouces;
  elemOld->computeExplicitSources(*geomCenter, *sourcesExplicit,
                                  numReadsExplicitSouces, 
                                  numWritesExplicitSouces
                                 );
  numReads  += numReadsExplicitSouces;
  numWrites += numWritesExplicitSouces;

  int numReadsImplicitSources, numWritesImplicitSources;
  elemOld->computeImplicitSources(*geomCenter, *sourcesImplicitOld,
			                            elemOld->tau,
                                  numReadsImplicitSources,
                                  numWritesImplicitSources
                                 );
  numReads  += numReadsImplicitSources;
  numWrites += numWritesImplicitSources;

  af::timer emhdGradientTimer = af::timer::start();
  if (params::viscosity || params::conduction)
  {
    int numReadsEMHDGradients, numWritesEMHDGradients;
    double dX[3];
    dX[0] = XCoords->dX1;
    dX[1] = XCoords->dX2;
    dX[2] = XCoords->dX3;
    elemOld->computeEMHDGradients(*geomCenter, dX,
                                  numReadsEMHDGradients,
                                  numWritesEMHDGradients
                                 );
    numReads  += numReadsEMHDGradients;
    numWrites += numWritesEMHDGradients;
  }
  double emhdGradientTime = af::timer::stop(emhdGradientTimer);
  
  af::timer divFluxTimer = af::timer::start();
  int numReadsDivFluxes, numWritesDivFluxes;
  computeDivOfFluxes(*primOld, numReadsDivFluxes, numWritesDivFluxes);
  numReads  += numReadsDivFluxes;
  numWrites += numWritesDivFluxes;
  double divFluxTime = af::timer::stop(divFluxTimer);

  /* Set a guess for prim */
  for (int var=0; var < vars::numFluidVars; var++)
  {
    prim->vars[var] = primOld->vars[var];
    numReads  += 1;
    numWrites += 1;
  }

  af::timer inductionEqnTimer = af::timer::start();
  cons->vars[vars::B1] = 
    consOld->vars[vars::B1] - 0.5*dt*divFluxes->vars[vars::B1];
  cons->vars[vars::B2] = 
    consOld->vars[vars::B2] - 0.5*dt*divFluxes->vars[vars::B2];
  cons->vars[vars::B3] = 
    consOld->vars[vars::B3] - 0.5*dt*divFluxes->vars[vars::B3];

  prim->vars[vars::B1] = cons->vars[vars::B1]/geomCenter->g;
  prim->vars[vars::B1].eval();
  prim->vars[vars::B2] = cons->vars[vars::B2]/geomCenter->g;
  prim->vars[vars::B2].eval();
  prim->vars[vars::B3] = cons->vars[vars::B3]/geomCenter->g;
  prim->vars[vars::B3].eval();
  double inductionEqnTime = af::timer::stop(inductionEqnTimer);

  primGuessPlusEps->vars[vars::B1] = prim->vars[vars::B1];
  primGuessPlusEps->vars[vars::B2] = prim->vars[vars::B2];
  primGuessPlusEps->vars[vars::B3] = prim->vars[vars::B3];

  primGuessLineSearchTrial->vars[vars::B1] = prim->vars[vars::B1];
  primGuessLineSearchTrial->vars[vars::B2] = prim->vars[vars::B2];
  primGuessLineSearchTrial->vars[vars::B3] = prim->vars[vars::B3];

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2 */
  jacobianAssemblyTime = 0.;
  lineSearchTime       = 0.;
  linearSolverTime     = 0.;
  af::timer solverTimer = af::timer::start();
  solve(*prim);
  double solverTime = af::timer::stop(solverTimer);

  /* Copy solution to primHalfStepGhosted. WARNING: Right now
   * primHalfStep->vars[var] points to prim->vars[var]. Might need to do a deep
   * copy. */
  for (int var=0; var < prim->numVars; var++)
  {
    primHalfStep->vars[var] = prim->vars[var];
    numReads  += 1;
    numWrites += 1;
  }
  af::timer halfStepCommTimer = af::timer::start();
  primHalfStep->communicate();
  double halfStepCommTime = af::timer::stop(halfStepCommTimer);

  af::timer halfStepDiagTimer = af::timer::start();
  halfStepDiagnostics(numReads,numWrites);
  double halfStepDiagTime = af::timer::stop(halfStepDiagTimer);
  /* Half step complete */

  double halfStepTime = af::timer::stop(halfStepTimer);

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "    ---Performance report--- \n");
  PetscPrintf(PETSC_COMM_WORLD, "     Boundary Conditions : %g secs, %g %\n",
                                 boundaryTime, boundaryTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Setting elemOld     : %g secs, %g %\n",
                                 elemOldTime, 
                                 elemOldTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Conserved vars Old  : %g secs, %g %\n",
                                 consOldTime, consOldTime/halfStepTime * 100
             );
  if (params::viscosity || params::conduction)
  {
    PetscPrintf(PETSC_COMM_WORLD, "     EMHD gradients      : %g secs, %g %\n",
                                 emhdGradientTime, 
                                 emhdGradientTime/halfStepTime * 100
               );
  }
  PetscPrintf(PETSC_COMM_WORLD, "     Divergence of fluxes: %g secs, %g %\n",
                                 divFluxTime, divFluxTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Nonlinear solver    : %g secs, %g %\n",
                                 solverTime, solverTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     |- Jacobian assembly: %g secs, %g %\n",
                                 jacobianAssemblyTime,
                                 jacobianAssemblyTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     |- CPU linear solver: %g secs, %g %\n",
                                 linearSolverTime, 
                                 linearSolverTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     |- Linesearch       : %g secs, %g %\n",
                                 lineSearchTime,
                                 lineSearchTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Induction equation  : %g secs, %g %\n",
                                 inductionEqnTime,
                                 inductionEqnTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Communication       : %g secs, %g %\n",
                                 halfStepCommTime,
                                 halfStepCommTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Diagnostics         : %g secs, %g %\n",
                                 halfStepDiagTime,
                                 halfStepDiagTime/halfStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Half step time      : %g secs\n\n",
                                 halfStepTime
             );

  /* Now take the full step */
  PetscPrintf(PETSC_COMM_WORLD, "  ---Full step--- \n");
  af::timer fullStepTimer = af::timer::start();

  currentStep = timeStepperSwitches::FULL_STEP;
  /* apply boundary conditions on primHalfStep */
  boundaryTimer = af::timer::start();
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *primHalfStep
                                     );
  setProblemSpecificBCs(numReads,numWrites);
  boundaryTime = af::timer::stop(boundaryTimer);

  af::timer elemHalfStepTimer = af::timer::start();
  elemHalfStep->set(*primHalfStep, *geomCenter,
                    numReadsElemSet, numWritesElemSet
                   );
  numReads  += numReadsElemSet;
  numWrites += numWritesElemSet; 
  double elemHalfStepTime = af::timer::stop(elemHalfStepTimer);

  elemHalfStep->computeExplicitSources(*geomCenter, *sourcesExplicit,
                                       numReadsExplicitSouces,
                                       numWritesExplicitSouces
                                      );
  numReads  += numReadsExplicitSouces;
  numWrites += numWritesExplicitSouces;

  elemOld->computeImplicitSources(*geomCenter, *sourcesImplicitOld,
			                            elemHalfStep->tau,
                                  numReadsImplicitSources,
                                  numWritesImplicitSources
                                 );
  numReads  += numReadsImplicitSources;
  numWrites += numWritesImplicitSources;

  emhdGradientTimer = af::timer::start();
  if (params::viscosity || params::conduction)
  {
    int numReadsEMHDGradients, numWritesEMHDGradients;
    double dX[3];
    dX[0] = XCoords->dX1;
    dX[1] = XCoords->dX2;
    dX[2] = XCoords->dX3;
    elemHalfStep->computeEMHDGradients(*geomCenter, dX,
                                       numReadsEMHDGradients,
                                       numWritesEMHDGradients
                                      );
    numReads  += numReadsEMHDGradients;
    numWrites += numWritesEMHDGradients;
  }
  emhdGradientTime = af::timer::stop(emhdGradientTimer);

  divFluxTimer = af::timer::start();
  computeDivOfFluxes(*primHalfStep, numReadsDivFluxes, numWritesDivFluxes);
  numReads  += numReadsDivFluxes;
  numWrites += numWritesDivFluxes;
  divFluxTime = af::timer::stop(divFluxTimer);

  inductionEqnTimer = af::timer::start();
  cons->vars[vars::B1] = 
    consOld->vars[vars::B1] - dt*divFluxes->vars[vars::B1];
  cons->vars[vars::B2] = 
    consOld->vars[vars::B2] - dt*divFluxes->vars[vars::B2];
  cons->vars[vars::B3] = 
    consOld->vars[vars::B3] - dt*divFluxes->vars[vars::B3];

  prim->vars[vars::B1] = cons->vars[vars::B1]/geomCenter->g;
  prim->vars[vars::B1].eval();
  prim->vars[vars::B2] = cons->vars[vars::B2]/geomCenter->g;
  prim->vars[vars::B2].eval();
  prim->vars[vars::B3] = cons->vars[vars::B3]/geomCenter->g;
  prim->vars[vars::B3].eval();
  inductionEqnTime = af::timer::stop(inductionEqnTimer);

  primGuessPlusEps->vars[vars::B1] = prim->vars[vars::B1];
  primGuessPlusEps->vars[vars::B2] = prim->vars[vars::B2];
  primGuessPlusEps->vars[vars::B3] = prim->vars[vars::B3];

  primGuessLineSearchTrial->vars[vars::B1] = prim->vars[vars::B1];
  primGuessLineSearchTrial->vars[vars::B2] = prim->vars[vars::B2];
  primGuessLineSearchTrial->vars[vars::B3] = prim->vars[vars::B3];

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2. NOTE: prim already has
   * primHalfStep as a guess */
  jacobianAssemblyTime = 0.;
  lineSearchTime       = 0.;
  linearSolverTime     = 0.;
  solverTimer = af::timer::start();
  solve(*prim);
  solverTime = af::timer::stop(solverTimer);

  /* Copy solution to primOldGhosted */
  for (int var=0; var < prim->numVars; var++)
  {
    primOld->vars[var] = prim->vars[var];
    numReads  += 1;
    numWrites += 1;
  }
  /* Compute diagnostics */
  af::timer fullStepCommTimer = af::timer::start();
  primOld->communicate();
  double fullStepCommTime = af::timer::stop(fullStepCommTimer);

  time += dt;
  af::timer fullStepDiagTimer = af::timer::start();
  fullStepDiagnostics(numReads,numWrites);
  double fullStepDiagTime = af::timer::stop(fullStepDiagTimer);

  /* done */
  double fullStepTime = af::timer::stop(fullStepTimer);
  double timeStepTime = af::timer::stop(timeStepTimer);

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "    ---Performance report--- \n");
  PetscPrintf(PETSC_COMM_WORLD, "     Boundary Conditions : %g secs, %g %\n",
                                 boundaryTime, boundaryTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Setting elemHalfStep: %g secs, %g %\n",
                                 elemHalfStepTime, 
                                 elemHalfStepTime/fullStepTime * 100
             );
  if (params::viscosity || params::conduction)
  {
    PetscPrintf(PETSC_COMM_WORLD, "     EMHD gradients      : %g secs, %g %\n",
                                 emhdGradientTime, 
                                 emhdGradientTime/fullStepTime * 100
               );
  }
  PetscPrintf(PETSC_COMM_WORLD, "     Divergence of fluxes: %g secs, %g %\n",
                                 divFluxTime, divFluxTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Nonlinear solver    : %g secs, %g %\n",
                                 solverTime, solverTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     |- Jacobian assembly: %g secs, %g %\n",
                                 jacobianAssemblyTime,
                                 jacobianAssemblyTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     |- CPU linear solver: %g secs, %g %\n",
                                 linearSolverTime, 
                                 linearSolverTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     |- Linesearch       : %g secs, %g %\n",
                                 lineSearchTime,
                                 lineSearchTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Induction equation  : %g secs, %g %\n",
                                 inductionEqnTime,
                                 inductionEqnTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Communication       : %g secs, %g %\n",
                                 fullStepCommTime,
                                 fullStepCommTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Diagnostics         : %g secs, %g %\n",
                                 fullStepDiagTime,
                                 fullStepDiagTime/fullStepTime * 100
             );
  PetscPrintf(PETSC_COMM_WORLD, "     Full step time      : %g secs\n\n",
                                 fullStepTime
             );
  PetscPrintf(PETSC_COMM_WORLD, "   ---Performance / proc : %g Zone cycles/sec/proc\n",
                                 prim->N1Local
                               * prim->N2Local
                               * prim->N3Local / timeStepTime
             );
  PetscPrintf(PETSC_COMM_WORLD, "   ---Total Performance  : %g Zone cycles/sec\n",
                                 N1 * N2 * N3 
                               / timeStepTime
             );
}

double timeStepper::computeDt(int &numReads, int &numWrites)
{
  // Time step control
  array minSpeedTemp,maxSpeedTemp;
  array minSpeed,maxSpeed;
  elemOld->computeMinMaxCharSpeeds(*geomCenter,
                                   directions::X1,
                                   minSpeedTemp, maxSpeedTemp,
                                   numReads,numWrites
                                  );
  minSpeedTemp = minSpeedTemp/XCoords->dX1;
  maxSpeedTemp = maxSpeedTemp/XCoords->dX1;
  maxSpeed     = af::max(maxSpeedTemp,af::abs(minSpeedTemp));

  if(params::dim>1)
  {
    elemOld->computeMinMaxCharSpeeds(*geomCenter,
                                     directions::X2,
                                     minSpeedTemp, maxSpeedTemp,
                                     numReads,numWrites
                                    );
    minSpeedTemp = minSpeedTemp/XCoords->dX2;
    maxSpeedTemp = maxSpeedTemp/XCoords->dX2;
    maxSpeed    += af::max(maxSpeedTemp,af::abs(minSpeedTemp));
  }

  if(params::dim>2)
  {
    elemOld->computeMinMaxCharSpeeds(*geomCenter,
                                     directions::X3,
                                     minSpeedTemp, maxSpeedTemp,
                                     numReads,numWrites);
    minSpeedTemp = minSpeedTemp/XCoords->dX3;
    maxSpeedTemp = maxSpeedTemp/XCoords->dX3;
    maxSpeed    += af::max(maxSpeedTemp,af::abs(minSpeedTemp));
  }
  array maxInvDt_af = af::max(af::max(af::max(maxSpeed,2),1),0);
  double maxInvDt = maxInvDt_af.host<double>()[0];

  /* Use MPI to find minimum over all processors */
  if (world_rank == 0) 
  {
    double temp; 
    for(int i=1;i<world_size;i++)
	  {
	    MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	    if( maxInvDt < temp)
      {
	      maxInvDt = temp;
      }
	  }
  }
  else
  {
    MPI_Send(&maxInvDt, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&maxInvDt,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  dt = params::CourantFactor/maxInvDt;
}
