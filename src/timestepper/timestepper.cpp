#include "timestepper.hpp"

timeStepper::timeStepper()
{
  geom                  = new geometry(0);
  geomGhosted           = new geometry(params::numGhost);

  prim                  = new grid(vars::dof, 0);                /* n+1   */
  primGhosted           = new grid(vars::dof, params::numGhost); /* n+1   */
  primHalfStepGhosted   = new grid(vars::dof, params::numGhost); /* n+1/2 */
  primOldGhosted        = new grid(vars::dof, params::numGhost); /* n     */

  cons                  = new grid(vars::dof, 0); /* n+1 */
  consOld               = new grid(vars::dof, 0); /* n   */

  sources               = new grid(vars::dof, 0); /* n     */
  sourcesHalfStep       = new grid(vars::dof, 0); /* n+1/2 */
  sourcesOld            = new grid(vars::dof, 0); /* n     */

  /* Fluxes at n or n+1/2 time step. Depends on context used */
  fluxesX1Ghosted       = new grid(vars::dof, params::numGhost); 
  fluxesX2Ghosted       = new grid(vars::dof, params::numGhost); 
  fluxesX3Ghosted       = new grid(vars::dof, params::numGhost); 

  divFluxes               = new grid(vars::dof, 0);

  elem                  = new fluidElement(*prim, geom, locations::CENTER); /* n+1 */
  elemOldGhosted        = new fluidElement(*primOldGhosted, geomGhosted,    /* n   */
                                           locations::CENTER
                                          );
  elemHalfStepGhosted   = new fluidElement(*primHalfStepGhosted, geomGhosted, 
                                           locations::CENTER                /* n+1/2*/
                                          );

  riemann = new riemannSolver(geomGhosted);

  nonLinSolver = new nonLinearSolver(computeResidual);

  SNES snes;
  if (   params::timeStepper==timeStepping::EXPLICIT
      || params::timeStepper==timeStepping::IMEX
     )
  {
    SNESSetDM(snes, prim->dm);
  }
  else if (params::timeStepper==timeStepping::IMPLICIT)
  {
    SNESSetDM(snes, prim->dmGhost);
  }

//  SNESSetFunction(snes, residual->globalVec, computeResidual, this);
  
  public:
    void timeStep(double dt);
}

timeStepper::~timeStepper()
{
  SNESDestroy(&snes);

  delete geom, geomGhosted;
  delete prim, primGhosted, primHalfStepGhosted, primOldGhosted;
  delete cons, consOld;
  delete sources, sourcesHalfStep, sourcesOld;
  delete residual;
  delete fluxesX1Ghosted,    fluxesX2Ghosted,    fluxesX3Ghosted;
  delete fluxesX1OldGhosted, fluxesX2OldGhosted, fluxesX3OldGhosted;
  delete fluxesX1HalfStepGhosted;
  delete fluxesX2HalfStepGhosted;
  delete fluxesX3HalfStepGhosted;
  delete elem, elemOldGhosted, elemHalfStepGhosted;
  delete riemann;
  delete nonLinSolver;
}





//PetscErrorCode computeResidual(SNES snes,
//                               Vec primVec,
//                               Vec residualVec,
//                               void *ptr
//                              )
//{
//  timeStepper *ts = (class timeStepper*)ptr;
//
//}


