#include "nonlinearsolver.hpp"

nonLinearSolver::nonLinearSolver(void (*computeResidual)(const grid &primGuess, 
                                                         grid &residual,
                                                         void *dataPtr
                                                        )
                                )
{
  this->computeResidual = computeResidual;

  primGuessLineSearchTrial  = new grid(vars::dof);
  primGuessPlusEps          = new grid(vars::dof);
  
  residual                  = new grid(vars::dof);
  residualPlusEps           = new grid(vars::dof);

  /* The grid data structure arranges data in Struct of Arrays format. Need to
   * rearrange to Array of Structs format in order to solve the linear system Ax
   * = b */
  residualSoA               = af::constant(0.,
                                           residual->vars[0].dims(0),
                                           residual->vars[0].dims(1),
                                           residual->vars[0].dims(2),
                                           vars::dof, f64
                                          );

  /* Jacobian \partial residual/ \prim in Struct of Arrays format */
  jacobianSoA               = af::constant(0, 
                                           residual->vars[0].dims(0),
                                           residual->vars[0].dims(1),
                                           residual->vars[0].dims(2),
                                           vars::dof * vars::dof, f64
                                          );

  /* RHS of Ax = b in Struct of Arrays format */
  bSoA                      = af::constant(0., 
                                           residual->vars[0].dims(0),
                                           residual->vars[0].dims(1),
                                           residual->vars[0].dims(2),
                                           vars::dof, f64
                                          );

  /* Correction dP_k in P_{k+1} = P_k + lambda*dP_k in Array of Structs format */
  deltaPrimAoS              = af::constant(0., 
                                           vars::dof,
                                           residual->vars[0].dims(0),
                                           residual->vars[0].dims(1),
                                           residual->vars[0].dims(2),
                                           f64
                                          );

  /* Steplength lambda in P_{k+1} = P_k + lambda*dP_k */ 
  stepLength                = af::constant(1., 
                                           residual->vars[0].dims(0),
                                           residual->vars[0].dims(1),
                                           residual->vars[0].dims(2),
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

