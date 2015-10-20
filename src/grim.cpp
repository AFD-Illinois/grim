#include "grim.h"
#include <yaml-cpp/yaml.h>
#include <arrayfire.h>

inline int DELTA(int const &mu, int const &nu)
{
  return (mu==nu ? 1 : 0);
}

using af::array;
using af::span;

const int NDIM = 4;
const int LOCATIONS = 7;
const int AXISYM_LOCATIONS = 5;

namespace boundaries
{
  enum
  {
    PERIODIC, OUTFLOW, MIRROR, DIRICHLET
  };
};

namespace locations
{
  enum
  {
    CENTER, LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
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
  int N1 = 128;
  int N2 = 128;
  int N3 = 128;
  int dim = 1;
  int numGhost = 2;

  int timeStepper = timeStepping::EXPLICIT;

  bool haveLocalIndicesBeenSet = 0;
  int N1Local, N2Local, N3Local;
  int iLocalStart, jLocalStart, kLocalStart;
  int iLocalEnd, jLocalEnd, kLocalEnd;

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
  int dof;
  int highOrderTermsConduction = 1;
  int highOrderTermsViscosity = 1;
  double adiabaticIndex = 4./3;

  int RHO = 0;
  int UU  = 1;
  int U1  = 2;
  int U2  = 3;
  int U3  = 4;
  int B1  = 5;
  int B2  = 6;
  int B3  = 7;
  int Q, DP;
};

class grid
{
  public:
    DM dm, dmGhost;
    Vec globalVec, localVec;

    int numVars, numGhost;

    DMBoundaryType boundaryLeft,  boundaryRight;
    DMBoundaryType boundaryTop,   boundaryBottom;
    DMBoundaryType boundaryFront, boundaryBack;

    int iLocalStart, iLocalEnd;
    int jLocalStart, jLocalEnd;
    int kLocalStart, kLocalEnd;

    array *vars;

    grid(int numVars, int numGhost);
    ~grid();

    void setVarsWithVec(Vec vec, bool needGhostZones);
};

class timeStepper
{
  public:
    SNES snes;
    grid *prim, *primOld;
    grid *cons, *consOld;
    grid *residual;

    fluidElement *elem, *elemOld;

    timeStepper();
    ~timeStepper();

    void timeStep();
};

void timeStepper::timeStep()
{
  elemOld->set(primOld->vars, geom, locations::CENTER);

  elemOld->computeFluxes(geom, 0, consOld->vars);
  elemOld->computeSources(geom, sourcesOld->vars);
  
  reconstruct(primOld->vars, 
              primLeft->vars, primRight->vars,
              directions::X1
             );

  riemannSolver(primRight->vars, primLeft->vars,
                geom, directions::X1, locations::LEFT,
                fluxesX1->vars
               );

}

PetscErrorCode computeResidual(SNES snes,
                               Vec primVec,
                               Vec residualVec,
                               void *ptr
                              )
{
  timeStepper *ts = (class timeStepper*)ptr;

}

timeStepper::timeStepper()
{

  SNES snes;
  prim = new grid(params::dof, params::numGhost);

  if (   params::timeStepper==timeStepping::EXPLICIT
      || params::timeStepper==timeStepping::IMEX
     )
  {
    SNESSetDM(snes, prim->dm);
  }
  else if (params::timeStepper==timeStepping::IMPLICIT)
  {
    SNESSetDM(snes, prim->dmGhost);
  }

  residual = new grid(params::dof, 0);
  primOld  = new grid(params::dof, params::numGhost);

  elem = new fluidElement(prim->vars, geom, locations::CENTER);
  elemOld = new fluidElement(primOld->vars, geom, locations::CENTER);

  SNESSetFunction(snes, residual->globalVec, computeResidual, this);

}

timeStepper::~timeStepper()
{
  SNESDestroy(&snes);
  
  delete elem;
  delete elemOld;

  delete prim;
  delete primOld;
  delete residual;
}

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
    params::Q   = 8;
    params::dof = 9;
  }
  else if (params::conduction==0 && params::viscosity==1)
  {
    params::DP  = 8;
    params::dof = 9;
  }
  else if (params::conduction==1 && params::viscosity==1)
  {
    params::Q   = 8;
    params::DP  = 9;
    params::dof = 10;
  }
  else if (params::conduction==0 && params::viscosity==0)
  {
    params::dof = 8;
  }


  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    timeStepper ts;
  }


  PetscFinalize();  
  return(0);
}

grid::grid(int numVars, int numGhost)
{
  this->numVars = numVars;
  this->numGhost = numGhost;

  /* Implementations for MIRROR, OUTFLOW in boundary.cpp and DIRICHLET in
   * problem.cpp */
  boundaryLeft  = DM_BOUNDARY_GHOSTED; boundaryRight  = DM_BOUNDARY_GHOSTED;
  boundaryTop   = DM_BOUNDARY_GHOSTED; boundaryBottom = DM_BOUNDARY_GHOSTED;
  boundaryFront = DM_BOUNDARY_GHOSTED; boundaryBack   = DM_BOUNDARY_GHOSTED;

  if (   params::boundaryLeft  == boundaries::PERIODIC 
      || params::boundaryRight == boundaries::PERIODIC
     )
  {
    boundaryLeft  = DM_BOUNDARY_PERIODIC;
    boundaryRight = DM_BOUNDARY_PERIODIC;
  }

  if (   params::boundaryTop    == boundaries::PERIODIC 
      || params::boundaryBottom == boundaries::PERIODIC
     )
  {
    boundaryTop    = DM_BOUNDARY_PERIODIC;
    boundaryBottom = DM_BOUNDARY_PERIODIC;
  }

  if (   params::boundaryFront == boundaries::PERIODIC 
      || params::boundaryBack  == boundaries::PERIODIC
     )
  {
    boundaryFront = DM_BOUNDARY_PERIODIC;
    boundaryBack  = DM_BOUNDARY_PERIODIC;
  }

  switch (params::dim)
  {
    case 1:
      DMDACreate1d(PETSC_COMM_WORLD, boundaryLeft, 
                   params::N1, numVars, 0, NULL,
                   &dm
                  );

      DMDACreate1d(PETSC_COMM_WORLD, boundaryLeft, 
                   params::N1, numVars, numGhost, NULL,
                   &dmGhost
                  );

      break;
  
    case 2:
      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, 0,
                   PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost,
                   PETSC_NULL, PETSC_NULL,
                   &dmGhost
                  );
      break;

    case 3:
      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2, params::N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, 0, 
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2, params::N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost, 
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dmGhost
                  );
      break;
  }

  DMCreateGlobalVector(dm, &globalVec);
  DMCreateLocalVector(dmGhost, &localVec);

  if (!params::haveLocalIndicesBeenSet)
  {
    DMDAGetCorners
      (dm,
       &params::iLocalStart, &params::jLocalStart, &params::kLocalStart,
       &params::N1Local,     &params::N2Local,     &params::N3Local
      );

    params::iLocalEnd = params::iLocalStart + params::N1Local;
    params::jLocalEnd = params::jLocalStart + params::N2Local;
    params::kLocalEnd = params::kLocalStart + params::N3Local;
  }

  vars = new array[numVars];

  if (numGhost > 0)
  {
    DMGlobalToLocalBegin(dmGhost, globalVec, INSERT_VALUES, localVec);
    DMGlobalToLocalEnd(dmGhost, globalVec, INSERT_VALUES, localVec);

    double *pointerToLocalVec;
    VecGetArray(localVec, &pointerToLocalVec);

    array varsCopiedFromVec = 
      array(numVars,
            params::N1Local + 2*numGhost, 
            params::N2Local + 2*numGhost,
            params::N3Local + 2*numGhost,
            pointerToLocalVec
           );

    for (int var=0; var<numVars; var++)
    {
      vars[var] = varsCopiedFromVec(var, span, span, span);
    }

    VecRestoreArray(localVec, &pointerToLocalVec);
  }
  else
  {
    double *pointerToGlobalVec;
    VecGetArray(globalVec, &pointerToGlobalVec);

    array varsCopiedFromVec = 
      array(numVars,
            params::N1Local, params::N2Local, params::N3Local,
            pointerToGlobalVec
           );

    for (int var=0; var<numVars; var++)
    {
      vars[var] = varsCopiedFromVec(var, span, span, span);
    }

    VecRestoreArray(globalVec, &pointerToGlobalVec);
  }
}

void grid::setVarsWithVec(Vec vec, bool needGhostZones)
{
  if (needGhostZones)
  {
    DMGlobalToLocalBegin(dmGhost, vec, INSERT_VALUES, localVec);
    DMGlobalToLocalEnd(dmGhost, vec, INSERT_VALUES, localVec);

    double *pointerToLocalVec;
    VecGetArray(localVec, &pointerToLocalVec);

    array varsCopiedFromVec = 
      array(numVars,
            params::N1Local + 2*params::numGhost, 
            params::N2Local + 2*params::numGhost,
            params::N3Local + 2*params::numGhost,
            pointerToLocalVec
           );

    for (int var=0; var<numVars; var++)
    {
      vars[var] = varsCopiedFromVec(var, span, span, span);
    }

    VecRestoreArray(globalVec, &pointerToLocalVec);
  }
  else
  {
    double *pointerToGlobalVec;
    VecGetArray(vec, &pointerToGlobalVec);

    array varsCopiedFromVec = 
      array(numVars,
            params::N1Local, params::N2Local, params::N3Local,
            pointerToGlobalVec
           );

    for (int var=0; var<numVars; var++)
    {
      vars[var] = varsCopiedFromVec(var, span, span, span);
    }

    VecRestoreArray(vec, &pointerToGlobalVec);
  }
}

grid::~grid()
{

  delete[] vars;

  VecDestroy(&globalVec);
  VecDestroy(&localVec);

  DMDestroy(&dm);
  DMDestroy(&dmGhost);
}

class geometry
{
  public:
    array xCoords[LOCATIONS][NDIM];
    array XCoords[LOCATIONS][NDIM];

    array alpha[AXISYM_LOCATIONS];
    array g[AXISYM_LOCATIONS];
    array gCov[AXISYM_LOCATIONS][NDIM][NDIM];
    array gCon[AXISYM_LOCATIONS][NDIM][NDIM];

    /* Connection coefficients only needed at CENTER */
    array gammaUpDownDown[NDIM][NDIM][NDIM];
};

class fluidElement
{
  private:
    const std::vector<std::vector<int>> indicesToLoopOver 
      = {{0},    {0, 1}, {0, 1}, {0, 2}, {0, 2}, {0, 3}, {0, 3}};
    /*   CENTER, LEFT,   RIGHT,  TOP,    BOTTOM, FRONT,  BACK*/
  public:
    int loc;

    /* fluidElement parameters */
    array tau, chi, nu;
    
    array rho, u, u1, u2, u3, B1, B2, B3;
    array pressure, temperature;
    array qTilde, deltaPTilde;
    array q, deltaP;
  
    array gammaLorentzFactor, uCon[NDIM], uCov[NDIM];
    array bSqr, bCon[NDIM], bCov[NDIM];
    
    array NUp[NDIM];
    array TUpDown[NDIM][NDIM];

    fluidElement();
    set(const array primVars[params::dof], 
        const geometry &geom, 
        const int location
       );
    void setFluidElementParameters(const geometry &geom);
    void computeFluxes(const geometry &geom, 
                       const int direction,
                       array fluxes[params::dof]
                      );
};

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = af::constant(1., params::N1Local, params::N2Local, params::N3Local);
  chi = af::constant(1., params::N1Local, params::N2Local, params::N3Local);
  nu  = af::constant(1., params::N1Local, params::N2Local, params::N3Local);
}

fluidElement::set(const array primVars[params::dof],
                  const geometry &geom,
                  const int location
                 )
{
  loc = location;

  rho = primVars[0] + params::rhoFloorInFluidElement;
  u   = primVars[1] + params::uFloorInFluidElement;
  u1  = primVars[2];
  u2  = primVars[3];
  u3  = primVars[4];
  B1  = primVars[5];
  B2  = primVars[6];
  B3  = primVars[7];

  setFluidElementParameters(geom);
  
  pressure    = (params::adiabaticIndex - 1.)*u;
  temperature = pressure/rho + params::temperatureFloorInFluidElement;

  if (params::conduction==1 && params::viscosity==0)
  {
    qTilde = primVars[8];

  }
  else if (params::conduction==0 && params::viscosity==1)
  {
    deltaPTilde = primVars[8];
  }
  else if (params::conduction==1 && params::viscosity==1)
  {
    qTilde      = primVars[8];
    deltaPTilde = primVars[9];
  }

  if (params::conduction==1)
  {
    if (params::highOrderTermsConduction==1)
    {
      q = qTilde * temperature * af::sqrt(rho*chi/tau);
    }
    else
    {
      q = qTilde;
    }
  }

  if (params::viscosity==1)
  {
    if (params::highOrderTermsViscosity == 1)
    {
      deltaP = deltaPTilde * af::sqrt(temperature * rho * nu / tau);
    }
    else
    {
      deltaP = deltaPTilde;
    }
  }

  gammaLorentzFactor = 
      sqrt(1 + geom.gCov[loc][1][1] * u1 * u1
             + geom.gCov[loc][2][2] * u2 * u2
             + geom.gCov[loc][3][3] * u3 * u3

           + 2*(  geom.gCov[loc][1][2] * u1 * u2
                + geom.gCov[loc][1][3] * u1 * u3
                + geom.gCov[loc][2][3] * u2 * u3
               )
          );
  
  uCon[0] = gamma/geom.alpha[loc];
  uCon[1] = u1 - gamma*geom.gCon[loc][0][1]*geom.alpha[loc];
  uCon[2] = u2 - gamma*geom.gCon[loc][0][2]*geom.alpha[loc];
  uCon[3] = u3 - gamma*geom.gCon[loc][0][3]*geom.alpha[loc];

  for (int mu=0; mu < NDIM; mu++)
  {
    uCov[mu] =  geom.gCov[loc][mu][0] * uCon[0]
              + geom.gCov[loc][mu][1] * uCon[1]
              + geom.gCov[loc][mu][2] * uCon[2]
              + geom.gCov[loc][mu][3] * uCon[3];
  }

  bCon[0] =  B1*uCov[1] + B2*uCov[2] + B3*uCov[3];

  bCon[1] = (B1 + bCon[0] * uCon[1])/uCon[0];
  bCon[2] = (B2 + bCon[0] * uCon[2])/uCon[0];
  bCon[3] = (B3 + bCon[0] * uCon[3])/uCon[0];

  for (int mu=0; mu < NDIM; mu++)
  {
    bCov[mu] =  geom.gCov[loc][mu][0] * bCon[0]
              + geom.gCov[loc][mu][1] * bCon[1]
              + geom.gCov[loc][mu][2] * bCon[2]
              + geom.gCov[loc][mu][3] * bCon[3];
  }

  bSqr =  bCon[0]*bCov[0] + bCon[1]*bCov[1]
        + bCon[2]*bCov[2] + bCon[3]*bCov[3] + params::bSqrFloorInFluidElement;


  for (int mu : indicesToLoopOver[loc])
  {
    NUp[mu] = rho * uCon[mu];

    for (int nu=0; nu < NDIM; nu++)
    {
      TUpDown[mu][nu] =   (rho + u + pressure + bSqr)*uCon[mu]*uCov[nu]
                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)
                        - bCon[mu] * bCov[nu]
                        + q/sqrt(bSqr) * (uCon[mu]*bCov[nu] + bCon[mu]*uCov[nu])
                        - deltaP       
                        * (  bCon[mu] * bCov[nu]/bSqr
                           + (1./3.)*(DELTA(mu, nu) + uCon[mu]*uCov[nu])
                          );
    }
  }

}

fluidElement::fluidElement(){}

void fluidElement::computeFluxes(const geometry &geom, 
                                 const int dir,
                                 array fluxes[params::dof]
                                )
{
  array g = geom.g[loc];

  fluxes[params::RHO] = g * NUp[dir];

  fluxes[params::UU]  = g * TUpDown[dir][0] + fluxes[params::RHO];

  fluxes[params::U1]  = g * TUpDown[dir][1];

  fluxes[params::U2]  = g * TUpDown[dir][2];

  fluxes[params::U3]  = g * TUpDown[dir][3];

  fluxes[params::B1]  = g * (bCon[1] * uCon[dir] - bCon[dir] * uCon[1]);

  fluxes[params::B2]  = g * (bCon[2] * uCon[dir] - bCon[dir] * uCon[2]);

  fluxes[params::B3]  = g * (bCon[3] * uCon[dir] - bCon[dir] * uCon[3]);

  if (params::conduction)
  {
    fluxes[params::Q] = g * (uCon[dir] * qTilde);
  }

  if (params::viscosity)
  {
    fluxes[params::DP] = g * (uCon[dir] * deltaPTilde);
  }
}

class riemannSolver
{
  fluidElement *elemLeft, *elemRight;

  array *fluxesLeft, *fluxesRight;
  array *consLeft,   *consRight;

  riemannSolver();
  ~riemannSolver();

  void solve(const array primVarsLeft[params::dof],
             const array primVarsRight[params::dof],
             const geometry &geom,
             const int dir,
             array fluxes[params::dof]
            );
};

riemannSolver::riemannSolver()
{
  elemLeft = new fluidElement;
}

void riemannSolver(const array primVarsLeft[params::dof],
                   const array primVarsRight[params::dof],
                   const geometry &geom,
                   const int dir,
                   array fluxes[params::dof]
                  )
{
  fluidElement elemLeft(primVarsLeft, geom, location);
  fluidElement elemRight(primVarsRight, geom, location);

  array fluxesLeft[params::dof], fluxesRight[params::dof];
  array consLeft[params::dof], consRight[params::dof];

  elemLeft.computeFluxes(geom, dir, fluxesLeft);
  elemLeft.computeFluxes(geom, 0,   consLeft);

  elemRight.computeFluxes(geom, dir, fluxesRight);
  elemRight.computeFluxes(geom, 0,   consRight);

  double cLaxFriedrichs = 1.;

  for (int var=0; var<params::dof; var++)
  {
    fluxes[var] = 0.5*(  fluxesLeft[var] + fluxesRight[var]
                       - cLaxFriedrichs * ( consRight[var] - consLeft[var])
                      );
              
  }
}
                                 
