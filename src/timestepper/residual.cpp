#include "timestepper.hpp"

void timeStepper::computeResidual(const grid &prim,
                                  grid &residual, 
                                  const bool computeExplicitTerms,
                                  int &numReads,
                                  int &numWrites
                                 )
{
  numReads = 0; numWrites = 0;
  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;
  elem->set(prim.vars, *geomCenter, numReadsElemSet, numWritesElemSet);
  elem->computeFluxes(*geomCenter, 0, cons->vars, 
                      numReadsComputeFluxes, numWritesComputeFluxes
                     );
  numReads  += numReadsElemSet  + numReadsComputeFluxes;
  numWrites += numWritesElemSet + numWritesComputeFluxes;

  if (currentStep == timeStepperSwitches::HALF_STEP)
  {
    if (computeExplicitTerms)
    {
      int numReadsExplicitSouces, numWritesExplicitSouces;
      elemOld->computeExplicitSources(*geomCenter, sourcesExplicit->vars,
                                      numReadsExplicitSouces, 
                                      numWritesExplicitSouces
                                     );
      numReads  += numReadsExplicitSouces;
      numWrites += numWritesExplicitSouces;
    }
    else
    {
      for (int var=0; var<vars::dof; var++)
      {
	      sourcesExplicit->vars[var] = 0.;
      }
    }

    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elemOld->computeImplicitSources(*geomCenter, sourcesImplicitOld->vars,
                                    numReadsImplicitSources,
                                    numWritesImplicitSources
                                   );
    elem->computeImplicitSources(*geomCenter, sourcesImplicit->vars,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources
                                );
    elemOld->computeTimeDerivSources(*geomCenter,
				                             *elemOld, *elem,
		                         		     dt/2,
                        				     sourcesTimeDer->vars,
                                     numReadsTimeDerivSources,
                                     numWritesTimeDerivSources
                                    );
    numReads  += 2*numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += 2*numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var<residual.numVars; var++)
    {
      residual.vars[var] = 
        (cons->vars[var] - consOld->vars[var])/(dt/2.)
  	  + divFluxes->vars[var]
      + sourcesExplicit->vars[var]
  	  + 0.5*(sourcesImplicitOld->vars[var] + sourcesImplicit->vars[var])
	    + sourcesTimeDer->vars[var];
    }
    /* Reads:
     * -----
     *  cons[var], consOld[var], divFluxes[var]       : 3*numVars
     *  sourcesExplicit[var], sourcesTimeDer[var]     : 2*numVars
     *  sourcesImplicitOld[var], sourcesImplicit[var] : 2*numVars
     *
     * Writes:
     * ------
     * residual[var] : numVars */
    numReads  += 7*residual.numVars;

    /* Normalization of the residual */
    if (params::conduction)
    {
	    if(params::highOrderTermsConduction)
      {
	      residual.vars[vars::Q] *=
          elemOld->temperature 
        * af::sqrt(elemOld->rho*elemOld->chi_emhd*elemOld->tau);

        numReads += 4;
      }
	    else
      {
	      residual.vars[vars::Q] *= elemOld->tau;
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
	    if(params::highOrderTermsViscosity)
      {
	      residual.vars[vars::DP] *=
          af::sqrt(   elemOld->rho*elemOld->nu_emhd
                    * elemOld->temperature*elemOld->tau
                  );
        numReads += 4;
      }
	    else
      {
	      residual.vars[vars::DP] *= elemOld->tau;
        numReads += 1;
      }
    }

  } /* End of timeStepperSwitches::HALF_STEP */

  else if (currentStep == timeStepperSwitches::FULL_STEP)
  {
    if (computeExplicitTerms)
    {
      int numReadsExplicitSouces, numWritesExplicitSouces;
      elemHalfStep->computeExplicitSources(*geomCenter, sourcesExplicit->vars,
                                           numReadsExplicitSouces,
                                           numWritesExplicitSouces
                                          );
      numReads  += numReadsExplicitSouces;
      numWrites += numWritesExplicitSouces;
    }
    else
    {
      for (int var=0; var<vars::dof; var++)
      {
	      sourcesExplicit->vars[var]=0.;
      }
    }

    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elemOld->computeImplicitSources(*geomCenter, sourcesImplicitOld->vars,
                                    numReadsImplicitSources,
                                    numWritesImplicitSources
                                   );
    elem->computeImplicitSources(*geomCenter, sourcesImplicit->vars,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources
                                );
    elemHalfStep->computeTimeDerivSources(*geomCenter,
	                                			  *elemOld, *elem,
					                                params::dt,
					                                sourcesTimeDer->vars,
                                          numReadsTimeDerivSources,
                                          numWritesTimeDerivSources
                                         );
    numReads  += 2*numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += 2*numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var<vars::dof; var++)
    {
      residual.vars[var] = 
        (cons->vars[var] - consOld->vars[var])/params::dt
    	+ divFluxes->vars[var]
	    + sourcesExplicit->vars[var]
	    + 0.5*(sourcesImplicitOld->vars[var] + sourcesImplicit->vars[var])
	    + sourcesTimeDer->vars[var];
    }
    /* Reads:
     * -----
     *  cons[var], consOld[var], divFluxes[var]       : 3*numVars
     *  sourcesExplicit[var], sourcesTimeDer[var]     : 2*numVars
     *  sourcesImplicitOld[var], sourcesImplicit[var] : 2*numVars
     *
     * Writes:
     * ------
     * residual[var] : numVars */
    numReads  += 7*residual.numVars;

    /* Normalization of the residual */
    if (params::conduction)
    {
	    if(params::highOrderTermsConduction)
      {
	      residual.vars[vars::Q] *= 
          elemHalfStep->temperature
        * af::sqrt(elemHalfStep->rho*elemHalfStep->chi_emhd*elemHalfStep->tau);

        numReads += 4;
      }
    	else
      {
	      residual.vars[vars::Q] *= elemHalfStep->tau;
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
      if (params::highOrderTermsViscosity)
      {
	      residual.vars[vars::DP] *= 
          af::sqrt(   elemHalfStep->rho*elemHalfStep->nu_emhd
                    * elemHalfStep->temperature*elemHalfStep->tau
                  );

        numReads += 4;
      }
	    else
      {
	      residual.vars[vars::DP] *= elemHalfStep->tau;
        numReads += 1;
      }
    }

  } /* End of timeStepperSwitches::FULL_STEP */

  for (int var=0; var < residual.numVars; var++)
  {
    residual.vars[var].eval();
  }
  numWrites += residual.numVars;
}
