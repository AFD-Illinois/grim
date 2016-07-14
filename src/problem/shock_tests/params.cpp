#include "../../params.hpp"
#include <arrayfire.h>

namespace params
{
  int numDevices = 1;

  int N1 = 2048;
  int N2 = 512;
  int N3 = 1;

  int dim = 1;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double InitialDt = 1e-5;
  double Time = 0.;
  double CourantFactor = 0.2;
  double maxDtIncrement = 1.01;
  double finalTime = 0.5;
  double WriteDataEveryDt = 0.1;
  int metric = metrics::MINKOWSKI;
  int restart = 0;
  std::string restartFile = "restartFile.h5";

  double X1Start = -2., X1End = 2.;
  double X2Start = 0., X2End = 1.;
  double X3Start = 0., X3End = 1.;

  int boundaryLeft   = boundaries::OUTFLOW;
  int boundaryRight  = boundaries::OUTFLOW;

  int boundaryTop    = boundaries::PERIODIC;
  int boundaryBottom = boundaries::PERIODIC;

  int boundaryFront  = boundaries::PERIODIC;
  int boundaryBack   = boundaries::PERIODIC;

  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;

  int conduction = 0;
  int viscosity  = 0;
  int highOrderTermsConduction = 0;
  int highOrderTermsViscosity = 0;
  double ConductionAlpha = 1.;
  double ViscosityAlpha = 1.;

  double adiabaticIndex = 4./3;

  std::string shockTest = "collision";

  double slopeLimTheta = 1;
  int reconstruction = reconstructionOptions::MINMOD;
  int riemannSolver  = riemannSolvers::LOCAL_LAX_FRIEDRICH;

  int maxNonLinearIter = 3;
  int maxLineSearchIters = 3;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-20;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;

  //Unused params - do we need to define them?
  double hSlope = 0.3;
  double blackHoleSpin = 0.9375;
  int DerefineThetaHorizon = 1;
  int DoCylindrify = 1;
  double X1cyl = 0.;
  double X2cyl = 1./N2;
};

namespace vars
{
  int Q   = 5;
  int DP  = 5 + params::conduction;
  int numFluidVars = 5 + params::conduction + params::viscosity;

  int B1  = 5 + params::conduction + params::viscosity;
  int B2  = 6 + params::conduction + params::viscosity;
  int B3  = 7 + params::conduction + params::viscosity;
  int dof = 8 + params::conduction + params::viscosity;
};
