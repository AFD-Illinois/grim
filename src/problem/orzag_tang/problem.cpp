#include "../problem.hpp"

void timeStepper::initialConditions(int &numReads, int &numWrites)
{
  double X1Center = (params::X1Start + params::X1End)/2.;
  double X2Center = (params::X2Start + params::X2End)/2.;

  array r = sqrt( pow((XCoords->vars[directions::X1] - X1Center), 2.) +
                  pow((XCoords->vars[directions::X2] - X2Center), 2.)
                );

  primOld->vars[vars::RHO] = 25./(36.*M_PI);
  primOld->vars[vars::U]   = 5./(12.*M_PI*(params::adiabaticIndex - 1.));
  //primOld->vars[vars::RHO] = 1. + af::exp(-r*r/0.01);
  //primOld->vars[vars::U]  = 1./(params::adiabaticIndex - 1.);

  array v1 = -0.5*af::sin(2*M_PI*XCoords->vars[directions::X2]);
  array v2 =  0.5*af::sin(2*M_PI*XCoords->vars[directions::X1]);
  //array v1 = 0.7 * elem->one;
  //array v2 = 0.7 * elem->one;
  array v3 =  0.*primOld->vars[vars::U3];
  array gamma = 1./af::sqrt(1. - v1*v1 - v2*v2 - v3*v3);

  primOld->vars[vars::U1] = gamma*v1;
  primOld->vars[vars::U2] = gamma*v2;
  primOld->vars[vars::U3] = gamma*v3;

  primOld->vars[vars::B1] = 
    -1./sqrt(4*M_PI) * af::sin(2*M_PI*XCoords->vars[directions::X2]);
  primOld->vars[vars::B2] = 
    1./sqrt(4*M_PI) * af::sin(4*M_PI*XCoords->vars[directions::X1]);
  primOld->vars[vars::B3] = 0.;

}

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi_emhd = 0.01 * one;
  nu_emhd  = 0.01 * one;
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{

}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{

}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{

}

void timeStepper::applyProblemSpecificFluxFilter(int &numReads,int &numWrites)
{

}
