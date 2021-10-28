#include "timestepper.hpp"

void timeStepper::computeResidual(const grid &primGuess,
                                  grid &residualGuess, 
                                  int &numReads,
                                  int &numWrites
                                 )
{
  numReads = 0; numWrites = 0;
  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;
  elem->set(primGuess, *geomCenter, numReadsElemSet, numWritesElemSet);
  elem->computeFluxes(0, *cons, 
                      numReadsComputeFluxes, numWritesComputeFluxes
                     );
  numReads  += numReadsElemSet  + numReadsComputeFluxes;
  numWrites += numWritesElemSet + numWritesComputeFluxes;

  if (currentStep == timeStepperSwitches::HALF_STEP)
  {
    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elem->computeImplicitSources(*sourcesImplicit,
                                 elemOld->tau,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources
                                );
    elemOld->computeTimeDerivSources(*elemOld, *elem,
                                     dt/2,
                                     *sourcesTimeDer,
                                     numReadsTimeDerivSources,
                                     numWritesTimeDerivSources
                                    );
    numReads  += numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var<residualGuess.numVars; var++)
    {
      residualGuess.vars[var] = 
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
     * residualGuess[var] : numVars */
    numReads  += 7*residualGuess.numVars;

    /* TODO: Need to make computeResidual() physics agnostic: 
     * the following code to normalize the residual needs to put somewhere else */
    /* Normalization of the residualGuess */
    if (params::conduction)
    {
      if(params::highOrderTermsConduction)
      {
        residualGuess.vars[vars::Q] *=
           elemOld->temperature 
         * af::sqrt(elemOld->rho*elemOld->chi_emhd*elemOld->tau);

        numReads += 4;
      }
      else
      {
        residualGuess.vars[vars::Q] *= elemOld->tau;
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
      if(params::highOrderTermsViscosity)
      {
        residualGuess.vars[vars::DP] *=
          af::sqrt(   elemOld->rho*elemOld->nu_emhd
                    * elemOld->temperature*elemOld->tau
                  );
        numReads += 4;
      }
      else
      {
        residualGuess.vars[vars::DP] *= elemOld->tau;
        numReads += 1;
      }
    }

  } /* End of timeStepperSwitches::HALF_STEP */

  else if (currentStep == timeStepperSwitches::FULL_STEP)
  {
    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elem->computeImplicitSources(*sourcesImplicit,
                                 elemHalfStep->tau,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources
                                );
    elemHalfStep->computeTimeDerivSources(*elemOld, *elem,
                                          dt,
                                          *sourcesTimeDer,
                                          numReadsTimeDerivSources,
                                          numWritesTimeDerivSources
                                         );
    numReads  += numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var < residualGuess.numVars; var++)
    {
      residualGuess.vars[var] = 
        (cons->vars[var] - consOld->vars[var])/dt
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
     * residualGuess[var] : numVars */
    numReads  += 7*residualGuess.numVars;

    /* Normalization of the residualGuess */
    if (params::conduction)
    {
      if(params::highOrderTermsConduction)
      {
        residualGuess.vars[vars::Q] *= 
          elemHalfStep->temperature
        * af::sqrt(elemHalfStep->rho*elemHalfStep->chi_emhd*elemHalfStep->tau);

        numReads += 4;
      }
      else
      {
        residualGuess.vars[vars::Q] *= elemHalfStep->tau;
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
      if (params::highOrderTermsViscosity)
      {
        residualGuess.vars[vars::DP] *= 
          af::sqrt(   elemHalfStep->rho*elemHalfStep->nu_emhd
                    * elemHalfStep->temperature*elemHalfStep->tau
                  );

        numReads += 4;
      }
      else
      {
        residualGuess.vars[vars::DP] *= elemHalfStep->tau;
        numReads += 1;
      }
    }

  } /* End of timeStepperSwitches::FULL_STEP */

  //Zero the residual in global ghost zones
  for (int var=0; var<residualGuess.numVars; var++) 
  {
    residualGuess.vars[var] *= residualMask;
  }

  numWrites += residualGuess.numVars;
}
