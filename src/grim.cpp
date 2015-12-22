#include "grim.hpp"

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

  if      (params::conduction==1 && params::viscosity==0)
  {
    vars::Q   = 8;
    vars::dof = 9;
  }
  else if (params::conduction==0 && params::viscosity==1)
  {
    vars::DP  = 8;
    vars::dof = 9;
  }
  else if (params::conduction==1 && params::viscosity==1)
  {
    vars::Q   = 8;
    vars::DP  = 9;
    vars::dof = 10;
  }
  else if (params::conduction==0 && params::viscosity==0)
  {
    vars::dof = 8;
  }


  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    
//    geometry geom(0);
//    geometry geomGhosted(params::numGhost);
//    grid primOldGhosted(vars::dof, params::numGhost);
//    grid consOld(vars::dof, 0);
//    grid consNew(vars::dof, 0);
//    grid fluxesX1OldGhosted(vars::dof, params::numGhost);
//    grid fluxesX2OldGhosted(vars::dof, params::numGhost);
//    grid fluxesX3OldGhosted(vars::dof, params::numGhost);
//
//    riemannSolver riemann(geomGhosted);
//
//    /* Initial conditions */
//    primOldGhosted.vars[vars::RHO] = 
//      1. + af::sin(2*M_PI*geomGhosted.xCoords[locations::CENTER][1]);
//    primOldGhosted.vars[vars::U]   = 1.;
//    primOldGhosted.vars[vars::U1]  = 0.;
//    primOldGhosted.vars[vars::U2]  = 0.;
//    primOldGhosted.vars[vars::U3]  = 0.;
//    primOldGhosted.vars[vars::B1]  = 0.;
//    primOldGhosted.vars[vars::B2]  = 0.;
//    primOldGhosted.vars[vars::B3]  = 0.;
//
//    riemann.solve(primOldGhosted, geomGhosted, directions::X1,
//                  fluxesX1OldGhosted);
//    riemann.solve(primOldGhosted, geomGhosted, directions::X2,
//                  fluxesX2OldGhosted);
//    riemann.solve(primOldGhosted, geomGhosted, directions::X3,
//                  fluxesX3OldGhosted);


  }

  PetscFinalize();  
  return(0);
}
