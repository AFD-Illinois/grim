#ifndef GRIM_PARAMS_H_
#define GRIM_PARAMS_H_

const int NDIM = 4;
const int LOCATIONS = 7;
const int AXISYM_LOCATIONS = 5;

namespace vars
{
  int RHO = 0;
  int U   = 1;
  int U1  = 2;
  int U2  = 3;
  int U3  = 4;
  int B1  = 5;
  int B2  = 6;
  int B3  = 7;
  int Q, DP;
  int dof;
};

namespace locations
{
  enum
  {
    CENTER, LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
  };
};

namespace directions
{
  enum
  {
    X1, X2, X3
  };
};


namespace gridParams
{
  bool haveGridParamsBeenSet = 0;
  int N1Local, N2Local, N3Local;
  int iLocalStart, jLocalStart, kLocalStart;
  int iLocalEnd,   jLocalEnd,   kLocalEnd;

  double dX1, dX2, dX3;
};

namespace boundaries
{
  enum
  {
    PERIODIC, OUTFLOW, MIRROR, DIRICHLET
  };
};

namespace metrics
{
  enum
  {
    MINKOWSKI, MODIFIED_KERR_SCHILD
  };
};

namespace timeStepping
{
  enum
  {
    EXPLICIT, IMEX, IMPLICIT
  };
};


namespace params
{
  int N1 = 8;
  int N2 = 8;
  int N3 = 8;
  int dim = 3;
  int numGhost = 2;

  int timeStepper = timeStepping::EXPLICIT;
  int metric = metrics::MINKOWSKI;
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

  int conduction = 0;
  int viscosity  = 0;
  int highOrderTermsConduction = 1;
  int highOrderTermsViscosity = 1;
  double adiabaticIndex = 4./3;

  double slopeLimTheta = 2;

};

#endif /* GRIM_PARAMS_H_ */
