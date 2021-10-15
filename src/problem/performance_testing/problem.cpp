#include "../problem.hpp"

void fluidElement::setFluidElementParameters()
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
}

void timeStepper::initialConditions(int &numReads,int &numWrites)
{
  array &rho = primOld->vars[vars::RHO];
  array &u   = primOld->vars[vars::U];
  array &u1  = primOld->vars[vars::U1];
  array &u2  = primOld->vars[vars::U2];
  array &u3  = primOld->vars[vars::U3];
  array &q   = primOld->vars[vars::Q];
  array &dP  = primOld->vars[vars::DP];

  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);
  int numGhostX1 = primOld->numGhostX1;
  int numGhostX2 = primOld->numGhostX2;
  int numGhostX3 = primOld->numGhostX3;
  int N1Local    = primOld->N1Local;
  int N2Local    = primOld->N2Local;
  int N3Local    = primOld->N3Local;

  for (int k=numGhostX3; k < N3Local+numGhostX3; k++)
  {
    for (int j=numGhostX2; j < N2Local+numGhostX2; j++)
    {
      for (int i=numGhostX1; i < N1Local+numGhostX1; i++)
      {
        double randNum = 0.5;
        
        PetscRandomGetValue(randNumGen, &randNum);
        rho(i,j,k) = 1 + 0.001*(randNum-0.5); // Needs to be positive

        PetscRandomGetValue(randNumGen, &randNum);
        u(i,j,k) = 1 + 0.001*(randNum-0.5); // same as above

        PetscRandomGetValue(randNumGen, &randNum);
        u1(i,j,k) = 0.001*(randNum-0.5);

        PetscRandomGetValue(randNumGen, &randNum);
        u2(i,j,k) = 0.001*(randNum-0.5);

        PetscRandomGetValue(randNumGen, &randNum);
        u3(i,j,k) = 0.001*(randNum-0.5);

        PetscRandomGetValue(randNumGen, &randNum);
        q(i,j,k) = 0.001*(randNum-0.5);

        PetscRandomGetValue(randNumGen, &randNum);
        dP(i,j,k) = 0.001*(randNum-0.5);
      }
    }
  }
  PetscRandomDestroy(&randNumGen);
}

void timeStepper::halfStepDiagnostics(int &numReadsTmp,int &numWritesTmp)
{
  //solve(*primOld);

  array jacobianAoS = af::reorder(jacobianSoA, 3, 0, 1, 2);

  /* RHS of Ax = b in Array of Structs format */
  array bAoS = -af::reorder(residualSoA, 3, 0, 1, 2);

  /* Now solve Ax = b using direct inversion, where
   * A = Jacobian
   * x = deltaPrim
   * b = -residual  */
  batchLinearSolve(jacobianAoS, bAoS, deltaPrimAoS);

  int numReads, numWrites;
  computeDivOfFluxes(*primHalfStep, numReads, numWrites);
  
  elemOld->set(*primOld, *geomCenter,
               numReads, numWrites
              );
  elemOld->computeFluxes(0, *consOld, 
                         numReads, numWrites
                        );
  af::sync();

  /* ======= End of Warm up ======== */

  af::timer linearSolverTimer = af::timer::start();
  for (int i=0; i < params::nIters; i++)
  {
    array jacobianAoS = af::reorder(jacobianSoA, 3, 0, 1, 2);
    array bAoS = -af::reorder(residualSoA, 3, 0, 1, 2);
    batchLinearSolve(jacobianAoS, bAoS, deltaPrimAoS);
    af::sync();
  }
  double linearSolverTime = af::timer::stop(linearSolverTimer);
  
  double zoneCyclesPerSec =  params::nIters
                           * (params::N1*params::N2*params::N3) 
                           / (linearSolverTime * 2. * params::maxNonLinearIter); // factor of 2 to account for half step and full step
  PetscPrintf(PETSC_COMM_WORLD, "\nTime taken for linear solver = %g secs, Zone cycles/sec = %g\n",
              linearSolverTime, zoneCyclesPerSec
             ); 

  af::sync();
  af::timer divFluxTimer = af::timer::start();
  for (int i=0; i < params::nIters; i++)
  {
    computeDivOfFluxes(*primHalfStep, numReads, numWrites);
    af::sync();
  }
  double divFluxTime = af::timer::stop(divFluxTimer);
  
  zoneCyclesPerSec =  params::nIters
                    * (params::N1*params::N2*params::N3) 
                    / (divFluxTime * 2.); // factor of 2 to account for half step and full step
  PetscPrintf(PETSC_COMM_WORLD, "\nTime taken for divFlux = %g secs, Zone cycles/sec = %g\n",
              divFluxTime, zoneCyclesPerSec
             ); 

  af::sync();
  af::timer elemSetTimer = af::timer::start();
  for (int i=0; i < params::nIters; i++)
  {
    elemOld->set(*primOld, *geomCenter,
                 numReads, numWrites
                );
    elemOld->computeFluxes(0, *consOld, 
                           numReads, numWrites
                          );
    af::eval(consOld->vars[vars::RHO],
             consOld->vars[vars::U],
             consOld->vars[vars::U1],
             consOld->vars[vars::U2],
             consOld->vars[vars::U3],
             consOld->vars[vars::B1],
             consOld->vars[vars::B2],
             consOld->vars[vars::B3],
             consOld->vars[vars::Q],
             consOld->vars[vars::DP]
            );
    af::sync();
  }
  double elemSetTime = af::timer::stop(elemSetTimer);
  
  int numCallsInRiemannSolver = 2;
  int numCallsInTimeStep = 2;
  int numCallsinComputeResidual = 1;
  int numCalls = 2 * (numCallsInTimeStep + (3*numCallsInRiemannSolver)); // 3 : one for each direction
  zoneCyclesPerSec =  params::nIters
                    * (params::N1*params::N2*params::N3) 
                    / (elemSetTime * numCalls);
  PetscPrintf(PETSC_COMM_WORLD, "\nTime taken for elem->Set(), elem->ComputeFluxes(0) [conserved vars] = %g secs, Zone cycles/sec = %g\n",
              elemSetTime, zoneCyclesPerSec
             ); 

  /* TODO: Check computeResidual() */
  
  PetscEnd();
}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{

}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{

}

void timeStepper::applyProblemSpecificFluxFilter(int &numReads,int &numWrites)
{

}


int timeStepper::CheckWallClockTermination()
{
  return 0;
}
