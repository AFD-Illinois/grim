#include "grim.hpp"

namespace vars
{
  int Q = 0;
  int DP = 0;
  int dof = 8;
};

namespace gridParams
{
  /* Set all vals to zero. Correct values will be set once a grid instance has
   * been initialized */
  bool haveGridParamsBeenSet = 0;
  int N1Local     = 0, N2Local      = 0, N3Local     = 0;
  int iLocalStart = 0, jLocalStart  = 0, kLocalStart = 0;
  int iLocalEnd   = 0, jLocalEnd    = 0, kLocalEnd   = 0;

  double dX1 = 0, dX2 = 0, dX3 = 0;
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
int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscPrintf(PETSC_COMM_WORLD, 
              "#### System info ####\n"
             );
  if (rank==0)
  {
    af::info();
  };

  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    
    geometry geom(0);
    geometry geomGhosted(params::numGhost);
    grid primOldGhosted(vars::dof, params::numGhost);
    grid consOld(vars::dof, 0);
    grid consNew(vars::dof, 0);
    grid fluxesX1OldGhosted(vars::dof, params::numGhost);
    grid fluxesX2OldGhosted(vars::dof, params::numGhost);
    grid fluxesX3OldGhosted(vars::dof, params::numGhost);

    riemannSolver riemann(geomGhosted);

    /* Initial conditions */
    primOldGhosted.vars[vars::RHO] = 
      1. + af::sin(2*M_PI*geomGhosted.xCoords[locations::CENTER][1]);
    primOldGhosted.vars[vars::U]   = 1.;
    primOldGhosted.vars[vars::U1]  = 0.;
    primOldGhosted.vars[vars::U2]  = 0.;
    primOldGhosted.vars[vars::U3]  = 0.;
    primOldGhosted.vars[vars::B1]  = 0.;
    primOldGhosted.vars[vars::B2]  = 0.;
    primOldGhosted.vars[vars::B3]  = 0.;

    riemann.solve(primOldGhosted, geomGhosted, directions::X1,
                  fluxesX1OldGhosted);
    riemann.solve(primOldGhosted, geomGhosted, directions::X2,
                  fluxesX2OldGhosted);
    riemann.solve(primOldGhosted, geomGhosted, directions::X3,
                  fluxesX3OldGhosted);


  }

  PetscFinalize();  
  return(0);
}
