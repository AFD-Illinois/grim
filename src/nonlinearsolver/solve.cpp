#include "timestepper.hpp"

void nonLinearSolver::solve(grid &primGuess)
{
  for (int nonLineariter=0;
       nonLineariter < params::maxNonLinearIter; nonLinearIter++
      )
  {
    computeResidual(primGuess, residual)

    for (int var=0; var < vars::dof; var++)
    {
      /* Need residualSoA to compute norms */
      residualSoA(span, span, span, var) = residual.vars[var];

      /* Initialize primGuessPlusEps. Needed to numerically assemble the
       * Jacobian */
      primGuessPlusEps.vars[var]         = primGuess.vars[var];
    }

    /* Still yet to put in code to compute global norm over mpi procs */
    double globalL2Norm =  af::norm(af::flat(residualSoA));
    printf("Nonliner iter = %d, error = %.15f\n", nonLinearIter, globalL2Norm);

    /* Assemble the Jacobian in Struct of Arrays format where the physics
     * operations are all vectorized */
    for (int row=0; row < vars::dof; row++)
    {
      /* Recommended value of Jacobian differencing parameter to achieve fp64
       * machine precision */
      double epsilon = 4e-8;

      primGuessPlusEps.vars[row]  = (1. + epsilon)*primGuess.vars[row]; 

      computeResidual(primGuessPlusEps, residualPlusEps)

      for (int column=0; column < vars::dof; column++)
      {
        jacobianSoA(span, span, span, column + vars::dof*row)
          = (  residualPlusEps.vars[column] 
             - residual.vars[column]
            )
            /(epsilon*primGuess.vars[row]);
      }

      /* reset */
      primGuessPlusEps.vars[row]  = primGuess.vars[row]; 
    }
    /* Jacobian assembly complete */

    /* Solve the linear system Jacobian * deltaPrim = -residual for the
     * correction deltaPrim */

    /* First assemble b. Need to reorder from SoA -> AoS. */
    for (int var=0; var < vars::dof; var++)
    {
      bSoA(span, span, span, var) = -residual.vars[var];
    }

    array jacobianAoS = af::reorder(jacobianSoA, 3, 0, 1, 2);
    array bAoS        = af::reorder(bSoA, 3, 0, 1, 2);

    /* Now solve Ax = b using direct inversion, where
     * A = Jacobian
     * x = deltaPrim
     * b = -residual 
     *
     * Currently inverting locally by looping over individual zones. Need to
     * call the batch function magma_dgesv_batched() from the MAGMA library
     * for optimal use on NVIDIA cards */
    for (int k=0; k<gridParams::N3Local; k++)
    {
      for (int j=0; j<gridParams::N2Local; j++)
      {
        for (int i=0; i<gridParams::N1Local; i++)
        {
          array A = af::moddims(jacobianAoS(span, i, j, k), 
                                vars::dof, vars::dof
                               );

          deltaPrimAoS(span, i, j, k) = af::solve(A, bAoS(span, i, j, k));
        }
      }
    }

    /* Done with the solve. Now rearrange from AoS -> SoA */
    array deltaPrimSoA = af::reorder(deltaPrimAoS, 1, 2, 3, 0);

    /* Quartic backtracking :
     * funcToMinimize \equiv 0.5 * normsL2 * normsL2
     */
    array f0      = 0.5 * af::sum(af::pow(residualSoA, 2.), 3);
    array fPrime0 = -2.*f0;
    
    array residualAoS = af::reorder(residualSoA, 3, 0, 1, 2);

    /* Start with a full step */
    stepLength = 1.;
    for (int lineSearchIter=0; 
         lineSearchIter < params::maxLineSearchIters; lineSearchIter++
        )
    {
      /* 1) First take the full step */
      for (int var=0; var<vars::dof; var++)
      {
        primGuessTrial.vars[var] =  
          primGuess.vars[var] + stepLength*deltaPrimSoA(span, span, span, var);
      } 

      /* ...and then compute the norm */
      computeResidual(primGuessTrial, residual);
      for (int var=0; var<vars::dof; var++)
      {
        residualSoA(span, span, span, var) = residual.vars[var];
      }
      array f1 = 0.5 * af::sum(af::pow(residualSoA, 2.), 3);

      /* We have 3 pieces of information:
       * a) f(0)
       * b) f'(0) 
       * c) f(1) 
       */
    
      /* NOT COMPLETE!!! */
      double alpha    = 1e-4;
      array condition = f1 > f0 + alpha*fPrime0;
      array lamda     = -fPrime0/(2.*(f1 - f0 - fPrime0));
      stepLength      = 1. - condition + condition*lamda;

//      array conditionIndices = where(condition > 0);
//      if (conditionIndices.elements() == 0)
//      {
//        break;
//      }
    }

    /* stepLength has now been set */
    for (int var=0; var<vars::dof; var++)
    {
      primGuess.vars[var] = 
        primGuess.vars[var] + stepLength*deltaPrimSoA(span, span, span, var);
    }
  }

}
