#include "physics.hpp"

riemannSolver::riemannSolver(const geometry &geom)
{
  primLeft  = new grid(vars::dof);
  primRight = new grid(vars::dof);

  elemLeft  = new fluidElement(*primLeft, geom, locations::LEFT);
  elemRight = new fluidElement(*primRight, geom, locations::LEFT);

  fluxLeft  = new grid(vars::dof);
  fluxRight = new grid(vars::dof);

  consLeft  = new grid(vars::dof);
  consRight = new grid(vars::dof);
}

riemannSolver::~riemannSolver()
{
  delete primLeft, primRight;
  delete elemLeft, elemRight;
  delete fluxLeft, fluxRight;
  delete consLeft, consRight;
}

void riemannSolver::solve(const grid &prim,
                          const geometry &geom,
                          const int dir,
                          grid &flux
                         )
{
  int location;
  int shiftX1, shiftX2, shiftX3;
  int fluxDirection;
  switch (dir)
  {
    case directions::X1:
      location = locations::LEFT;
      fluxDirection = 1;
      shiftX1  = 1;
      shiftX2  = 0;
      shiftX3  = 0;
      break;

    case directions::X2:
      location = locations::BOTTOM;
      fluxDirection = 2;
      shiftX1  = 0;
      shiftX2  = 1;
      shiftX3  = 0;
      break;

    case directions::X3:
      location = locations::BACK;
      fluxDirection = 3;
      shiftX1  = 0;
      shiftX2  = 0;
      shiftX3  = 1;
      break;
  }

  reconstruction::reconstruct(prim, dir, *primLeft, *primRight);

  elemLeft->set(*primRight, geom, location);
  elemRight->set(*primLeft, geom, location);

  elemLeft->computeFluxes(geom, fluxDirection, *fluxLeft);
  elemLeft->computeFluxes(geom, 0,             *consLeft);

  elemRight->computeFluxes(geom, fluxDirection, *fluxRight);
  elemRight->computeFluxes(geom, 0,             *consRight);

  double cLaxFriedrichs = 1.;

  for (int var=0; var<vars::dof; var++)
  {
    flux.vars[var] = 0.5*(  af::shift(fluxLeft->vars[var], shiftX1, shiftX2, shiftX3) 
                          + fluxRight->vars[var]
                          - cLaxFriedrichs 
                          * (  consRight->vars[var]
			       - af::shift(consLeft->vars[var], shiftX1, shiftX2, shiftX3)
                            )
                         );
  }
            
}
