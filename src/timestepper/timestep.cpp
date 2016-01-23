#include "timestepper.hpp"

void timeStepper::computeDivOfFluxes(const grid &primFlux)
{
  switch (params::dim)
  {
    case 1:
      riemann->solve(primFlux, *geom, directions::X1, *fluxesX1);

      for (int var=0; var<vars::dof; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxesX1->dX1);

        array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);

        divFluxes->vars[var] = dFluxX1_dX1;
      }

      break;

    case 2:
      riemann->solve(primFlux, *geom, directions::X1, *fluxesX1);
      riemann->solve(primFlux, *geom, directions::X2, *fluxesX2);

      for (int var=0; var<vars::dof; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxesX1->dX1);
        array filterX2 = array(1, 3, 1, 1, filter1D)/(fluxesX2->dX2);

        array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);
        array dFluxX2_dX2 = convolve(fluxesX2->vars[var], filterX2);

        divFluxes->vars[var] = dFluxX1_dX1 + dFluxX2_dX2;
      }

      break;

    case 3:
      riemann->solve(primFlux, *geom, directions::X1, *fluxesX1);
      riemann->solve(primFlux, *geom, directions::X2, *fluxesX2);
      riemann->solve(primFlux, *geom, directions::X3, *fluxesX3);

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

      break;
  }
}

void timeStepper::timeStep()
{
  /* First take a half step */
  currentStep = timeStepperSwitches::HALF_STEP;

  grid primTile(params::tileSizeX1,
                params::tileSizeX2, 
                1+2*params::numGhost
               );

  grid XCoordsTile(params::tileSizeX1, 
                   params::tileSizeX2,
                   1+2*params::numGhost
                  );

  geometry geomTile(XCoordsTile);



  for (int k=0; k < params::N3; k++)
  {
    gfor (af::seq jTile, numTilesX2)
    {
      for (int iTile; iTile < numTilesX1; iTile++)
      {
        indices->setTile(iTile, jTile, k);
        setXCoords(indices, locations::CENTER, *XCoordsTile);

        primTile->copy(*prim, *XCoordsTile);
        geomTile->copy(*geom, *XCoordsTile);

        elemTile->set(*primTile, *geomTile);

        elemTile->computeFluxes(*geomTile, 0, *consOldTile);

      }
    }
  }

  setProblemSpecificBCs();

  elemOld->set(*primOld, *geom, locations::CENTER);
  elemOld->computeFluxes(*geom, 0, *consOld);
  if(params::viscosity || params::conduction)
    elemOld->computeEMHDGradients(*geom);
  
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
  halfStepDiagnostics();
  elemHalfStep->set(*primHalfStep, *geom, locations::CENTER);
  /* apply boundary conditions on primHalfStepGhosted */
  /* Half step complete */

  /* Now take the full step */
  currentStep = timeStepperSwitches::FULL_STEP;
  setProblemSpecificBCs();

  if(params::viscosity || params::conduction)
    elemHalfStep->computeEMHDGradients(*geom);
  computeDivOfFluxes(*primHalfStep);

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2. NOTE: prim already has
   * primHalfStep as a guess */
  solve(*prim);

  /* Copy solution to primOldGhosted */
  for (int var=0; var<vars::dof; var++)
  {
    primOld->vars[var] = prim->vars[var];
  }
  fullStepDiagnostics();

  /* Compute diagnostics */

  /* done */

}

void timeStepper::computeResidual(const grid &primGuess, grid &residualGuess, const bool ComputeExplicitTerms)
{
  elem->set(primGuess, *geom, locations::CENTER);
  elem->computeFluxes(*geom, 0, *cons);

  if (currentStep == timeStepperSwitches::HALF_STEP)
  {
    int useImplicitSources = 0;
    if(ComputeExplicitTerms)
      elemOld->computeExplicitSources(*geom, 
				      *sourcesE
				      );
    else
      for (int var=0; var<vars::dof; var++)
	sourcesE->vars[var]=0.;
    elemOld->computeImplicitSources(*geom,
				    *sourcesIOld);
    elem->computeImplicitSources(*geom,
				 *sourcesINew);
    elemOld->computeTimeDerivSources(*geom,
				     *elemOld,
				     *elem,
				     params::dt/2,
				     *sourcesDT);
    for (int var=0; var<vars::dof; var++)
    {
      residualGuess.vars[var] = (  cons->vars[var] 
				   - consOld->vars[var]
				   )/(params::dt/2.)
	+ divFluxes->vars[var]
	+ sourcesE->vars[var]
	+0.5*(sourcesIOld->vars[var]+sourcesINew->vars[var])
	+sourcesDT->vars[var];
    }
    //Normalization of the residual
    for (int var=0; var<vars::dof; var++)
      residualGuess.vars[var] = residualGuess.vars[var]/geom->g[locations::CENTER];
    if(params::conduction)
      {
	if(params::highOrderTermsConduction)
	  residualGuess.vars[vars::Q] = residualGuess.vars[vars::Q]*elemOld->temperature*af::sqrt(elemOld->rho*elemOld->chi_emhd*elemOld->tau);
	else
	  residualGuess.vars[vars::Q] = residualGuess.vars[vars::Q]*elemOld->tau;
      }
    if(params::viscosity)
      {
	if(params::highOrderTermsViscosity)
	  residualGuess.vars[vars::DP] = residualGuess.vars[vars::DP]*af::sqrt(elemOld->rho*elemOld->nu_emhd*elemOld->temperature*elemOld->tau);
	else
	  residualGuess.vars[vars::DP] = residualGuess.vars[vars::DP]*elemOld->tau;
      }
  }
  else if (currentStep == timeStepperSwitches::FULL_STEP)
  {
    int useImplicitSources = 0;
    if(ComputeExplicitTerms)
      elemHalfStep->computeExplicitSources(*geom, 
					   *sourcesE
					   );
    else
      for (int var=0; var<vars::dof; var++)
	sourcesE->vars[var]=0.;
    elemOld->computeImplicitSources(*geom,
				    *sourcesIOld);
    elem->computeImplicitSources(*geom,
				 *sourcesINew);
    elemHalfStep->computeTimeDerivSources(*geom,
					  *elemOld,
					  *elem,
					  params::dt,
					  *sourcesDT);
    for (int var=0; var<vars::dof; var++)
    {
      residualGuess.vars[var] = (  cons->vars[var] 
				   - consOld->vars[var]
				   )/params::dt
	+ divFluxes->vars[var]
	+ sourcesE->vars[var]
	+0.5*(sourcesIOld->vars[var]+sourcesINew->vars[var])
	+ sourcesDT->vars[var];
    }
   //Normalization of the residual
    for (int var=0; var<vars::dof; var++)
      residualGuess.vars[var] = residualGuess.vars[var]/geom->g[locations::CENTER];
    if(params::conduction)
      {
	if(params::highOrderTermsConduction)
	  residualGuess.vars[vars::Q] = residualGuess.vars[vars::Q]*elemHalfStep->temperature*af::sqrt(elemHalfStep->rho*elemHalfStep->chi_emhd*elemHalfStep->tau);
	else
	  residualGuess.vars[vars::Q] = residualGuess.vars[vars::Q]*elemHalfStep->tau;
      }
    if(params::viscosity)
      {
	if(params::highOrderTermsViscosity)
	  residualGuess.vars[vars::DP] = residualGuess.vars[vars::DP]*af::sqrt(elemHalfStep->rho*elemHalfStep->nu_emhd*elemHalfStep->temperature*elemHalfStep->tau);
	else
	  residualGuess.vars[vars::DP] = residualGuess.vars[vars::DP]*elemHalfStep->tau;
      }
  }

  //This is one way to guarantee that ghost zones are synced...
  residualGuess.communicate();
}
