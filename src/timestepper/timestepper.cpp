#include "timestepper.hpp"

timeStepper::timeStepper()
{
  geom          = new geometry();

  prim          = new grid(vars::dof);  /* n+1   */
  primHalfStep  = new grid(vars::dof);  /* n+1/2 */
  primOld       = new grid(vars::dof);  /* n     */

  cons          = new grid(vars::dof); /* n+1 */
  consOld       = new grid(vars::dof); /* n   */

  sources         = new grid(vars::dof); /* n     */
  sourcesHalfStep = new grid(vars::dof); /* n+1/2 */
  sourcesOld      = new grid(vars::dof); /* n     */

  /* Fluxes at n or n+1/2 time step. Depends on context used */
  fluxesX1  = new grid(vars::dof); 
  fluxesX2  = new grid(vars::dof); 
  fluxesX3  = new grid(vars::dof); 

  divFluxes = new grid(vars::dof);

  elem          = new fluidElement(*prim,         *geom, locations::CENTER); /* n+1   */
  elemOld       = new fluidElement(*primOld,      *geom, locations::CENTER); /* n     */
  elemHalfStep  = new fluidElement(*primHalfStep, *geom, locations::CENTER); /* n+1/2 */

  riemann = new riemannSolver(*geom);

  nonLinSolver = new nonLinearSolver(computeResidual);
  
}

timeStepper::~timeStepper()
{
  delete geom;
  delete prim, primHalfStep, primOld;
  delete cons, consOld;
  delete sources, sourcesHalfStep, sourcesOld;
  delete fluxesX1, fluxesX2, fluxesX3;
  delete elem, elemOld, elemHalfStep;
  delete riemann;
  delete nonLinSolver;
}
