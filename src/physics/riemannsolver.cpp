#include "physics.hpp"

riemannSolver::riemannSolver(const geometry &geom)
{
  primLeft  = new grid(vars::dof, params::numGhost);
  primRight = new grid(vars::dof, params::numGhost);

  elemLeft  = new fluidElement(*primLeft, geom, locations::LEFT);
  elemRight = new fluidElement(*primRight, geom, locations::LEFT);

  fluxLeft  = new grid(vars::dof, 0);
  fluxRight = new grid(vars::dof, 0);

  consLeft  = new grid(vars::dof, 0);
  consRight = new grid(vars::dof, 0);
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
  switch (dir)
  {
    case directions::X1:
      location = locations::LEFT;
      break;

    case directions::X2:
      location = locations::BOTTOM;
      break;

    case directions::X3:
      location = locations::BACK;
      break;
  }

  reconstruct(prim, dir, *primLeft, *primRight);

  elemLeft->set(*primRight, geom, location);
  elemRight->set(*primLeft, geom, location);

  elemLeft->computeFluxes(geom, dir, *fluxLeft);
  elemLeft->computeFluxes(geom, 0,   *consLeft);

  elemRight->computeFluxes(geom, dir, *fluxRight);
  elemRight->computeFluxes(geom, 0,   *consRight);

  double cLaxFriedrichs = 1.;

  for (int var=0; var<vars::dof; var++)
  {
    flux.vars[var] = 0.5*(  fluxLeft->vars[var] + fluxRight->vars[var]
                          - cLaxFriedrichs 
                          * (consRight->vars[var] - consLeft->vars[var])
                         );
  }
            
}

array minmod(array &x, array &y, array &z)
{
  array minOfAll = af::min(
                           af::min(af::abs(x), af::abs(y)), 
                           af::abs(z)
                          );

  return 0.25 * abs(sign(x) + sign(y) ) * (sign(x) + sign(z) ) * minOfAll;
}

void riemannSolver::reconstruct(const grid &prim,
                                const int dir,
                                grid &primLeft,
                                grid &primRight
                               )
{
  double filter1D[]  = {1,-1, 0, /* Forward difference */
                        0, 1,-1  /* Backward difference */
                       };
  array filter;
  switch (dir)
  {
    case directions::X1:
      filter =  array(3, 1, 1, 2, filter1D)/gridParams::dX1;
      break;

    case directions::X2:
      filter =  array(1, 3, 1, 2, filter1D)/gridParams::dX2;
      break;

    case directions::X3:
      filter =  array(1, 1, 3, 2, filter1D)/gridParams::dX3;
      break;
  }

  for(int var=0; var<vars::dof; var++)
  {
    array dvar_dX = convolve(prim.vars[var], filter);

    array forwardDiff  = dvar_dX(span, span, span, 0);
    array backwardDiff = dvar_dX(span, span, span, 1);
    array centralDiff  = backwardDiff + forwardDiff;

    array left   = params::slopeLimTheta * backwardDiff;
    array center = 0.5 * centralDiff;
    array right  = params::slopeLimTheta * forwardDiff;

    array slope  =  minmod(left, center, right);

    primLeft.vars[var]  = prim.vars[var] - 0.5*slope;
    primRight.vars[var] = prim.vars[var] + 0.5*slope;
  }

}
