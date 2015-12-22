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
  double dt = 0.01;
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

    for (int var=0; var<vars::dof; var++)
    {
      double filter1D[] = {1, -1, 0};
      
      array filterX1 = array(3, 1, 1, 1, filter1D)/gridParams::dX1;
      array filterX2 = array(1, 3, 1, 1, filter1D)/gridParams::dX2;
      array filterX3 = array(1, 1, 3, 1, filter1D)/gridParams::dX3;

      array dFluxX1_dX1 = convolve(fluxesX1OldGhosted.vars[var], filterX1);
      array dFluxX2_dX2 = convolve(fluxesX2OldGhosted.vars[var], filterX2);
      array dFluxX3_dX3 = convolve(fluxesX3OldGhosted.vars[var], filterX3);

      af::seq domainX1(params::numGhost, params::N1 + params::numGhost - 1);
      af::seq domainX2(params::numGhost, params::N2 + params::numGhost - 1);
      af::seq domainX3(params::numGhost, params::N3 + params::numGhost - 1);

      consNew.vars[var] =  consOld.vars[var] 
                         + params::dt*(  dFluxX1_dX1(domainX1, domainX2, domainX3)
                                       + dFluxX2_dX2(domainX1, domainX2, domainX3)
                                       + dFluxX3_dX3(domainX1, domainX2, domainX3)
                                      );;
    }

  }

  PetscFinalize();  
  return(0);
}
