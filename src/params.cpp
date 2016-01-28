#include "params.hpp"

namespace vars
{
  int Q = 8;
  int DP = 9;
  int dof = 10;
};

namespace params
{
  int N1 = 48;
  int N2 = 48;
  int N3 = 48;

  int dim = 3;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double dt = .5/N1;
  double Time = 0.;
  double finalTime = 0.5;
  int metric = metrics::MODIFIED_KERR_SCHILD;
  double hSlope = 0.3;
  double blackHoleSpin = 0.9375;

  double X1Start = 0., X1End = 1.;
  double X2Start = 0., X2End = 1.;
  double X3Start = 0., X3End = 1.;

  int boundaryLeft   = boundaries::PERIODIC;
  int boundaryRight  = boundaries::PERIODIC;

  int boundaryTop    = boundaries::PERIODIC;
  int boundaryBottom = boundaries::PERIODIC;

  int boundaryFront  = boundaries::PERIODIC;
  int boundaryBack   = boundaries::PERIODIC;

  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;

  int conduction = 1;
  int viscosity  = 1;
  int highOrderTermsConduction = 1;
  int highOrderTermsViscosity = 1;
  double adiabaticIndex = 4./3;

  double slopeLimTheta = 2;
  int reconstruction = reconstructionOptions::WENO5;

  int maxNonLinearIter = 10;
  int maxLineSearchIters = 10;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-10;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;
  
};
