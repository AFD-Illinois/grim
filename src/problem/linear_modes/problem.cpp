#include "../problem.hpp"

void fluidElement::setFluidElementParameters()
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
}

void timeStepper::initialConditions(int &numReads,int &numWrites)
{
  array xCoords[3];
  for(int d=0;d<3;d++)
    xCoords[d]=XCoords->vars[d];
  geomCenter->XCoordsToxCoords(XCoords->vars,xCoords);

  array cphi = af::cos(  params::k1*xCoords[directions::X1]
                       + params::k2*xCoords[directions::X2]
                      );

  array sphi = af::sin(  params::k1*xCoords[directions::X1]
                       + params::k2*xCoords[directions::X2]
                      );

  /* Initial conditions */

  //Alfven Wave
  /*primOld->vars[vars::RHO] = 1.;
  primOld->vars[vars::U]   = 2.;
  primOld->vars[vars::U1]  = 0.;
  primOld->vars[vars::U2]  = params::Aw*0.462905090215*cphi;
  primOld->vars[vars::U3]  = 0.;
  primOld->vars[vars::B1]  = 0.01;
  primOld->vars[vars::B2]  = params::Aw*0.886407850514*cphi;
  primOld->vars[vars::B3]  = 0.;*/

  //Sound wave
  /*primOld->vars[vars::RHO] = 1.+params::Aw*0.345991032308*cphi;
  primOld->vars[vars::U]   = 2.+params::Aw*0.922642752822*cphi;
  primOld->vars[vars::U1]  = 0.-params::Aw*0.170354208129*cphi; 
  primOld->vars[vars::U2]  = 0.; 
  primOld->vars[vars::U3]  = 0.; 
  primOld->vars[vars::B1]  = 0.01; 
  primOld->vars[vars::B2]  = 0.;
  primOld->vars[vars::B3]  = 0.;*/

  //Full EMHD mode (from grim2D)
  primOld->vars[vars::RHO] = 1.;
  primOld->vars[vars::U]   = 2.;
  primOld->vars[vars::U1]  = 0.;
  primOld->vars[vars::U2]  = 0.; 
  primOld->vars[vars::U3]  = 0.; 
  primOld->vars[vars::B1]  = 0.1; 
  primOld->vars[vars::B2]  = 0.3;
  primOld->vars[vars::B3]  = 0.;
  primOld->vars[vars::Q]   = 0.;
  primOld->vars[vars::DP]  = 0.;

  primOld->vars[vars::RHO] += params::Aw*cphi*(-0.518522524082246)
                    +params::Aw*sphi*0.1792647678001878;

  primOld->vars[vars::U]   += params::Aw*cphi*0.5516170736393813;

  primOld->vars[vars::U1]  += params::Aw*cphi*0.008463122479547856
                    +params::Aw*sphi*(-0.011862022608466367);

  primOld->vars[vars::U2]  += params::Aw*cphi*(-0.16175466371870734)
                    +params::Aw*sphi*(0.034828080823603294);

  primOld->vars[vars::B1]  += params::Aw*cphi*(-0.05973794979640743)
                    +params::Aw*sphi*0.03351707506150924;

  primOld->vars[vars::B2]  += params::Aw*cphi*0.02986897489820372
                    -params::Aw*sphi*0.016758537530754618;

  primOld->vars[vars::Q]   += params::Aw*cphi*0.5233486841539436
                    -params::Aw*sphi*0.04767672501939603;

  primOld->vars[vars::DP]  += params::Aw*cphi*0.2909106062057657
                    -params::Aw*sphi*0.02159452055336572;

  array pressure    = (params::adiabaticIndex - 1.)*primOld->vars[vars::U];
  array temperature = pressure/primOld->vars[vars::RHO];
  array soundSpeed  = af::sqrt(params::adiabaticIndex*pressure
                               / (primOld->vars[vars::RHO]
                                  + params::adiabaticIndex * primOld->vars[vars::U]
                                 )
                              );

  if (params::highOrderTermsConduction)
  {
    primOld->vars[vars::Q] = 
      primOld->vars[vars::Q] / temperature
    / af::sqrt(primOld->vars[vars::RHO] * params::ConductionAlpha * soundSpeed * soundSpeed);
  }

  if (params::highOrderTermsViscosity)
  {
    primOld->vars[vars::DP] = 
      primOld->vars[vars::DP]
    / af::sqrt(  temperature * primOld->vars[vars::RHO] * params::ViscosityAlpha
               * soundSpeed * soundSpeed
              );
  }

  for (int var=0; var < numVars; var++)
  {
    primOld->vars[var].eval();
  }
  af::sync();
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{

}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
  /* Compute the errors for the different modes */

  array xCoords[3];
  for(int d=0;d<3;d++)
    xCoords[d]=XCoords->vars[d];
  geomCenter->XCoordsToxCoords(XCoords->vars,xCoords);

  //Alfven wave
  /*cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
    +0.0328124176673*params::Time);
    array u2an = params::Aw*0.462905090215*cphi;
    double error = af::norm(af::flat((ts.primOld->vars[vars::U2]
    - u2an)));*/
  
  //MHD Sound wave
  /*cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
    +3.09362659024*params::Time);
    array rhoan = 1.+params::Aw*0.345991032308*cphi;
    double error = af::norm(af::flat((ts.primOld->vars[vars::RHO]
    - rhoan)));*/
  
  //EMHD Sound wave
  array cphi = af::cos(  params::k1*xCoords[directions::X1]
       + params::k2*xCoords[directions::X2]
       + params::Omega*time
       );
  
  array sphi = af::sin(  params::k1*xCoords[directions::X1]
       + params::k2*xCoords[directions::X2]
       + params::Omega*time
       );
  
  array rhoAnalytic = 1. + (  params::Aw*cphi*(-0.518522524082246)
            + params::Aw*sphi*0.1792647678001878
            )*exp(params::Gamma*time);
  
  array dPAnalytic = 0. + ( params::Aw*cphi*0.2909106062057657
          -params::Aw*sphi*0.02159452055336572
          )*exp(params::Gamma*time);

  array uAnalytic = 2. + params::Aw*cphi*0.5516170736393813*exp(params::Gamma*time);

  array u1Analytic = (  params::Aw*cphi*0.008463122479547856
                      + params::Aw*sphi*(-0.011862022608466367)
                     ) * exp(params::Gamma*time);

  array u2Analytic = (  params::Aw*cphi*(-0.16175466371870734)
                      + params::Aw*sphi*(0.034828080823603294)
                     ) * exp(params::Gamma*time);

  array qAnalytic =  (  params::Aw*cphi*0.5233486841539436
                      - params::Aw*sphi*0.04767672501939603
                     ) * exp(params::Gamma*time);

  array B1Analytic = 0.1 + (  params::Aw*cphi*(-0.05973794979640743)
                            + params::Aw*sphi*0.03351707506150924
                           ) * exp(params::Gamma*time);

  array B2Analytic = 0.3 + (  params::Aw*cphi*0.02986897489820372
                            - params::Aw*sphi*0.016758537530754618
                           ) * exp(params::Gamma*time);
  
  elemOld->set(*primOld, *geomCenter, numReads, numWrites);

  double errorRho = af::norm(af::flat(elemOld->rho     - rhoAnalytic));
  double errorU   = af::norm(af::flat(elemOld->u       - uAnalytic));
  double errorU1  = af::norm(af::flat(elemOld->u1      - u1Analytic));
  double errorU2  = af::norm(af::flat(elemOld->u2      - u2Analytic));
  double errorB1  = af::norm(af::flat(elemOld->B1      - B1Analytic));
  double errorB2  = af::norm(af::flat(elemOld->B2      - B2Analytic));
  double errorQ   = af::norm(af::flat(elemOld->q       - qAnalytic));
  double errordP  = af::norm(af::flat(elemOld->deltaP  - dPAnalytic));

  
  errorRho = errorRho/N1/N2/N3;
  errorU   = errorU/N1/N2/N3;
  errorU1  = errorU1/N1/N2/N3;
  errorU2  = errorU2/N1/N2/N3;
  errorQ   = errorQ/N1/N2/N3;
  errordP  = errordP/N1/N2/N3;
  errorB1  = errorB1/N1/N2/N3;
  errorB2  = errorB2/N1/N2/N3;

  PetscPrintf(PETSC_COMM_WORLD, "\n");
  PetscPrintf(PETSC_COMM_WORLD, "    ---Linear Modes Diagnostics---\n");
  PetscPrintf(PETSC_COMM_WORLD, "     Error in rho = %e\n", errorRho);
  PetscPrintf(PETSC_COMM_WORLD, "     Error in u   = %e\n", errorU  );
  PetscPrintf(PETSC_COMM_WORLD, "     Error in u1  = %e\n", errorU1 );
  PetscPrintf(PETSC_COMM_WORLD, "     Error in u2  = %e\n", errorU2 );
  PetscPrintf(PETSC_COMM_WORLD, "     Error in q   = %e\n", errorQ  );
  PetscPrintf(PETSC_COMM_WORLD, "     Error in dP  = %e\n", errordP );
  PetscPrintf(PETSC_COMM_WORLD, "     Error in B1  = %e\n", errorB1 );
  PetscPrintf(PETSC_COMM_WORLD, "     Error in B2  = %e\n", errorB2 );
}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{

}

void timeStepper::applyProblemSpecificFluxFilter(int &numReads,int &numWrites)
{

}


int timeStepper::CheckWallClockTermination()
{
  return 0;
}
