#include "timestepper.hpp"

void timeStepper::computeDivOfFluxes(const grid &prim)
{
  riemann->solve(prim, *geom, directions::X1, *fluxesX1);
  riemann->solve(prim, *geom, directions::X2, *fluxesX2);
  riemann->solve(prim, *geom, directions::X3, *fluxesX3);

  for (int var=0; var<vars::dof; var++)
  {
    double filter1D[] = {1, -1, 0}; /* Forward difference */
    
    array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxesX1->dX1);
    array filterX2 = array(1, 3, 1, 1, filter1D)/(fluxesX2->dX2);
    array filterX3 = array(1, 1, 3, 1, filter1D)/(fluxesX3->dX3);

    array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);
    array dFluxX2_dX2 = convolve(fluxesX2->vars[var], filterX2);
    array dFluxX3_dX3 = convolve(fluxesX3->vars[var], filterX3);

    divFluxes->vars[var] = dFluxX1_dX1 + dFluxX2_dX2 + dFluxX3_dX3;
  }

}

void timeStepper::timeStep(double dt)
{
  /* First take a half step */
  currentStep = timeStepperSwitches::HALF_STEP;

  elemOld->set(*primOld, *geom, locations::CENTER);
  elemOld->computeFluxes(*geom, 0, *consOld);

  computeDivOfFluxes(*primOld);

  /* Set a guess for prim */
  for (int var=0; var<vars::dof; var++)
  {
    prim->vars[var] = primOld->vars[var];
  }

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2 */
  solve(*prim);

  /* Copy solution to primHalfStepGhosted. WARNING: Right now
   * primHalfStep->vars[var] points to prim->vars[var]. Might need to do a deep
   * copy. */
  for (int var=0; var<vars::dof; var++)
  {
    primHalfStep->vars[var] = prim->vars[var];
  }

  /* apply boundary conditions on primHalfStepGhosted */
  /* Half step complete */

  /* Now take the full step */
  currentStep = timeStepperSwitches::FULL_STEP;

  computeDivOfFluxes(*primHalfStep);

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2. NOTE: prim already has
   * primHalfStep as a guess */
  solve(*prim);

  /* Copy solution to primOldGhosted */
  for (int var=0; var<vars::dof; var++)
  {
    primOld->vars[var] = prim->vars[var];
  }

  /* Compute diagnostics */

  /* done */

}

void timeStepper::computeResidual(const grid &prim, grid &residual)
{

  elem->set(prim, *geom, locations::CENTER);
  elem->computeFluxes(*geom, 0, *cons);

  if (currentStep == timeStepperSwitches::HALF_STEP)
  {
    int useImplicitSources = 0;
    elem->computeSources(*geom, 
                         *elemOld, *elemOld,
                         params::dt/2., useImplicitSources,
                         *sourcesOld
                        );

    for (int var=0; var<vars::dof; var++)
    {
      residual.vars[var] = (  cons->vars[var] 
                            - consOld->vars[var]
                           )/(params::dt/2.)
                          + divFluxes->vars[var]
                          + sourcesOld->vars[var];
    }

  }
  else if (currentStep == timeStepperSwitches::FULL_STEP)
  {
    int useImplicitSources = 0;
    elem->computeSources(*geom, 
                         *elemOld, *elemHalfStep, 
                         params::dt, useImplicitSources,
                         *sourcesOld
                        );

    for (int var=0; var<vars::dof; var++)
    {
      residual.vars[var] = (  cons->vars[var] 
                            - consOld->vars[var]
                           )/params::dt
                          + divFluxes->vars[var]
                          + sourcesOld->vars[var];
    }
  }


}
