#include "nonlinearsolver.hpp"

nonLinearSolver::nonLinearSolver(void (*computeResidual)(const grid &primGuess, 
                                                         grid &residual
                                                        )
                                )
{
  this->computeResidual = computeResidual;

  primGuessLineSearchTrial  = new grid(vars::dof, 0);
  primGuessPlusEps          = new grid(vars::dof, 0);
  
  residual                  = new grid(vars::dof, 0);
  residualPlusEps           = new grid(vars::dof, 0);

  /* The grid data structure arranges data in Struct of Arrays format. Need to
   * rearrange to Array of Structs format in order to solve the linear system Ax
   * = b */
  residualSoA               = af::constant(0.,
                                           gridParams::N1Local,
                                           gridParams::N2Local,
                                           gridParams::N3Local,
                                           vars::dof, f64
                                          );

  /* Jacobian \partial residual/ \prim in Struct of Arrays format */
  jacobianSoA               = af::constant(0, 
                                           gridParams::N1Local,
                                           gridParams::N2Local,
                                           gridParams::N3Local,
                                           vars::dof * vars::dof, f64
                                          );

  /* RHS of Ax = b in Struct of Arrays format */
  bSoA                      = af::constant(0., 
                                           gridParams::N1Local,
                                           gridParams::N2Local,
                                           gridParams::N3Local,
                                           vars::dof, f64
                                          );

  /* Correction dP_k in P_{k+1} = P_k + lambda*dP_k in Array of Structs format */
  deltaPrimAoS              = af::constant(0., 
                                           vars::dof,
                                           gridParams::N1Local,
                                           gridParams::N2Local,
                                           gridParams::N3Local,
                                           f64
                                          );

  /* Steplength lambda in P_{k+1} = P_k + lambda*dP_k */ 
  steplength                = af::constant(1., 
                                           gridParams::N1Local,
                                           gridParams::N2Local,
                                           gridParams::N3Local,
                                           f64
                                          );
  
  
}

nonLinearSolver::~nonLinearSolver()
{
  delete primGuessLineSearchTrial;
  delete primGuessPlusEps;
  
  delete residual;
  delete residualPlusEps;
}

