#include "timestepper.hpp"

timeStepper::timeStepper()
{
  geom          = new geometry();

  prim          = new grid(vars::dof);  /* n+1   */
  primHalfStep  = new grid(vars::dof);  /* n+1/2 */
  primOld       = new grid(vars::dof);  /* n     */

  cons          = new grid(vars::dof); /* n+1 */
  consOld       = new grid(vars::dof); /* n   */

  sourcesE         = new grid(vars::dof); /* Explicit     */
  sourcesIOld      = new grid(vars::dof); /* Implicit at n */
  sourcesINew      = new grid(vars::dof); /* Implicit at n+1 */
  sourcesDT        = new grid(vars::dof); /* TimeDeriv     */

  /* Fluxes at n or n+1/2 time step. Depends on context used */
  fluxesX1  = new grid(vars::dof); 
  fluxesX2  = new grid(vars::dof); 
  fluxesX3  = new grid(vars::dof); 

  divFluxes = new grid(vars::dof);

  elem          = new fluidElement(*prim,         *geom, locations::CENTER); /* n+1   */
  elemOld       = new fluidElement(*primOld,      *geom, locations::CENTER); /* n     */
  elemHalfStep  = new fluidElement(*primHalfStep, *geom, locations::CENTER); /* n+1/2 */

  riemann = new riemannSolver(*geom);

  /* Data structures needed for the nonlinear solver */
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

timeStepper::~timeStepper()
{
  delete geom;
  delete prim, primHalfStep, primOld;
  delete cons, consOld;
  delete sourcesE, sourcesINew, sourcesIOld, sourcesDT;
  delete fluxesX1, fluxesX2, fluxesX3;
  delete elem, elemOld, elemHalfStep;
  delete riemann;

  delete primGuessLineSearchTrial;
  delete primGuessPlusEps;
  delete residual;
  delete residualPlusEps;
}
