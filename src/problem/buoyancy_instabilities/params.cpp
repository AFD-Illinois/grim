#include "../../params.hpp"
#include <arrayfire.h>

namespace params
{
  int numDevices = 4;

  int N1 = 64;
  int N2 = 64;
  int N3 = 128;
  int dim = 2;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double InitialDt = 0.01;
  double CourantFactor = 0.5;
  double maxDtIncrement = 1.3;
  double Time = 0;
  double finalTime = 10.;

  int restart = 0;
  std::string restartFile = "restartFile.h5";

  int ObserveEveryNSteps = 100;
  int StepNumber = 0;

  double adiabaticIndex = 5./3;
  double blackHoleSpin = 0.;
  double InnerEdgeRadius = 6.;
  double PressureMaxRadius = 12.;
  double MinPlasmaBeta  = 15.;
  double MagneticLoops = 1;
  double Adiabat = 0.001;

  double Rin = 100*(1.+sqrt(1.-blackHoleSpin*blackHoleSpin));
  double Rout = 300.;
  
  // Grid parameters
  int metric = metrics::MODIFIED_KERR_SCHILD;
  double hSlope = 1.;
  int DerefineThetaHorizon = 1;
  int DoCylindrify = 1;
  double X1cyl = log(4.*Rin);
  double X2cyl = 3./N2;

  
  double X1Start = log(Rin), X1End = log(Rout);
  double X2Start = 0.+1.e-8, X2End = 1.-1.e-8;
  double X3Start = 0., X3End = M_PI*1.;

  int boundaryLeft   = boundaries::OUTFLOW;
  int boundaryRight  = boundaries::OUTFLOW;

  int boundaryTop    = boundaries::MIRROR;
  int boundaryBottom = boundaries::MIRROR;

  int boundaryFront  = boundaries::PERIODIC;
  int boundaryBack   = boundaries::PERIODIC;

  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;
  
  //Atmosphere parameters
  // Floors are Ampl*pow(radius,power)
  double RhoFloorAmpl = 1.e-3;
  double UFloorAmpl = 1.e-5;
  double RhoFloorSlope = -1.5;
  double UFloorSlope = -2.5;
  // Floors for magnetically dominated regions
  double BsqrOverRhoMax = 10.;
  double BsqrOverUMax = 500.;

  int conduction = 1;
  int viscosity  = 0;
  int highOrderTermsConduction = 1.;
  int highOrderTermsViscosity = 1.;
  double ConductionAlpha = 1.;
  double ViscosityAlpha = 1.;
  double ConductionClosureFactor = 1.;
  double ViscosityClosureFactor = 1.;

  double slopeLimTheta = 2;
  int reconstruction = reconstructionOptions::MINMOD;
  int riemannSolver  = riemannSolvers::LOCAL_LAX_FRIEDRICH;

  int maxNonLinearIter = 3;
  int maxLineSearchIters = 3;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-6;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;
  
  double InitialPerturbationAmplitude = 4e-2;
  double ObserveEveryDt = 1.;
  double WriteDataEveryDt = 20.;
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
