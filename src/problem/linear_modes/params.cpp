#include "../../params.hpp"
#include <arrayfire.h>

namespace params
{
  int numDevices = 2;

  int N1 = 64;
  int N2 = 64;
  int N3 = 64;

  int dim = 3;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double InitialDt = .5/N1;
  double maxDtIncrement = 1.3;
  double Time = 0.;
  double CourantFactor = 0.9;
  double finalTime = 0.1;
  double MaxWallTime = 3600*23.5;
  int metric = metrics::MINKOWSKI;
  int restart = 0;
  std::string restartFile = "restartFile.h5";
  std::string restartFileName = "restartFileName.txt";
  std::string restartFileTime = "restartFileTime.txt";

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
  double ConductionAlpha = 1.;
  double ViscosityAlpha = 1.;

  double adiabaticIndex = 4./3;
  double Aw = 1.e-8;
  double k1 = 2.*M_PI;
  double k2 = 4.*M_PI;
  double Gamma = - 0.5533585207638141;
  double Omega = - 3.6262571286888425;

  double slopeLimTheta = 2;
  int reconstruction = reconstructionOptions::MINMOD;
  int riemannSolver  = riemannSolvers::HLL;

  int maxNonLinearIter = 3;
  int maxLineSearchIters = 3;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-20;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;

  // Linear solver options
  int linearSolver = linearSolvers::GPU_BATCH_SOLVER;

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

namespace dumpVars
{
  int Q   = 5;
  int DP  = 5 + params::conduction;
  int numFluidVars = 5 + params::conduction + params::viscosity;

  int B1  = 5 + params::conduction + params::viscosity;
  int B2  = 6 + params::conduction + params::viscosity;
  int B3  = 7 + params::conduction + params::viscosity;
  int dof = 8 + params::conduction + params::viscosity;
};
