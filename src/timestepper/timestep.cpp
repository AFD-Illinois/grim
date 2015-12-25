#include "timestepper.hpp"

void timeStepper::computeDivOfFluxes(const grid &primGhosted)
{
  riemann->solve(primOldGhosted, geomGhosted, directions::X1,
                 fluxesX1Ghosted
                );
  riemann->solve(primOldGhosted, geomGhosted, directions::X2,
                 fluxesX2Ghosted
                );
  riemann->solve(primOldGhosted, geomGhosted, directions::X3,
                 fluxesX3Ghosted
                );

  for (int var=0; var<vars::dof; var++)
  {
    double filter1D[] = {1, -1, 0}; /* Forward difference */
    
    array filterX1 = array(3, 1, 1, 1, filter1D)/gridParams::dX1;
    array filterX2 = array(1, 3, 1, 1, filter1D)/gridParams::dX2;
    array filterX3 = array(1, 1, 3, 1, filter1D)/gridParams::dX3;

    array dFluxX1_dX1 = convolve(fluxesX1Ghosted.vars[var], filterX1);
    array dFluxX2_dX2 = convolve(fluxesX2Ghosted.vars[var], filterX2);
    array dFluxX3_dX3 = convolve(fluxesX3Ghosted.vars[var], filterX3);

    divFluxes.vars[var] = 
        dFluxX1_dX1(gridParams::domainX1, 
                    gridParams::domainX2, 
                    gridParams::domainX3
                   )
      + dFluxX2_dX2(gridParams::domainX1, 
                    gridParams::domainX2, 
                    gridParams::domainX3
                   )
      + dFluxX3_dX3(gridParams::domainX1, 
                    gridParams::domainX2, 
                    gridParams::domainX3
                   );
  }

}

void timeStepper::timeStep(double dt)
{
  /* First take a half step */
  currentStep = timeStepperSwitches::HALF_STEP;

  elemOldGhosted->set(primOldGhosted, geomGhosted, locations::CENTER);
  elemOldGhosted->computeFluxes(geomGhosted, 0, consOldGhosted);

  computeDivOfFluxes(primOldGhosted);

  /* Set a guess for prim */
  for (int var=0; var<vars::dof; var++)
  {
    prim.vars[var] = primOldGhosted.vars[var]
                      (gridParams::domainX1,
                       gridParams::domainX2,
                       gridParams::domainX3
                      );
  }

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2 */
  nonLinSolver->solve(prim);

  /* Copy solution to primHalfStepGhosted */
  for (int var=0; var<vars::dof; var++)
  {
    primHalfStepGhosted.vars[var]
      (gridParams::domainX1,
       gridParams::domainX2,
       gridParams::domainX3
      ) 
    = prim.vars[var];
  }

  /* apply boundary conditions on primHalfStepGhosted */
  /* Half step complete */

  /* Now take the full step */
  currentStep = timeStepperSwitches::FULL_STEP;

  computeDivOfFluxes(primHalfStepGhosted);

  /* Set a guess for prim */
  for (int var=0; var<vars::dof; var++)
  {
    prim.vars[var] = primHalfStepGhosted.vars[var]
                      (gridParams::domainX1,
                       gridParams::domainX2,
                       gridParams::domainX3
                      );
  }

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2 */
  nonLinSolver->solve(prim);

  /* Copy solution to primOldGhosted */
  for (int var=0; var<vars::dof; var++)
  {
    primOldGhosted.vars[var]
      (gridParams::domainX1,
       gridParams::domainX2,
       gridParams::domainX3
      ) 
    = prim.vars[var];
  }

  /* Compute diagnostics */

  /* done */

}

void computeResidual(const grid &prim, grid &residual, void *dataPtr)
{
  timeStepper *ts = static_cast<timeStepper *>dataPtr;

  ts->elem->set(prim, geom, locations::CENTER);
  ts->elem->computeFluxes(geom, 0, cons);

  if (ts->currentStep == timeStepperSwitches::HALF_STEP)
  {
    ts->elem->computeSources(geom, 
                             elemOldGhosted, elemOldGhosted, 
                             params::dt/2., sourcesOld
                            );

    for (int var=0; var<vars::dof; var++)
    {
      residual.vars[var] = (  ts->cons.vars[var] 
                            - ts->consOldGhosted.vars[var]
                                    (gridParams::domainX1,
                                     gridParams::domainX2,
                                     gridParams::domainX3
                                    );
                            )/(ts->dt/2.)
                          + ts->divFluxes.vars[var]
                          + ts->sourcesOld.vars[var];
    }

  }
  else if (ts->currentStep == timeStepperSwitches::FULL_STEP)
  {
    ts->elem->computeSources(geom, 
                             elemOldGhosted, elemHalfStepGhosted, 
                             params::dt, sourcesOld
                            );

    for (int var=0; var<vars::dof; var++)
    {
      residual.vars[var] = (  ts->cons.vars[var] 
                            - ts->consOldGhosted.vars[var]
                                    (gridParams::domainX1,
                                     gridParams::domainX2,
                                     gridParams::domainX3
                                    );
                            )/ts->dt
                          + ts->divFluxes.vars[var]
                          + ts->sourcesOld.vars[var];
    }



  }


}
