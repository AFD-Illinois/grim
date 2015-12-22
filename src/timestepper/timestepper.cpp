#include "timestepper.hpp"

timeStepper::timeStepper(const geometry &geom)
{

  prim        = new grid(vars::dof, params::numGhost);
  primOld     = new grid(vars::dof, params::numGhost);

  cons        = new grid(vars::dof, 0);
  consOld     = new grid(vars::dof, 0);

  sourcesOld  = new grid(vars::dof, 0);
  residual    = new grid(vars::dof, 0);

  switch (params::dim)
  {
    case 1:
      fluxesX1 = new grid(vars::dof, 1);
      break;

    case 2:
      fluxesX1 = new grid(vars::dof, 1);
      fluxesX2 = new grid(vars::dof, 1);
      break;

    case 3:
      fluxesX1 = new grid(vars::dof, 1);
      fluxesX2 = new grid(vars::dof, 1);
      fluxesX3 = new grid(vars::dof, 1);
      break;
  }

  elem    = new fluidElement(*prim, geom, locations::CENTER);
  elemOld = new fluidElement(*primOld, geom, locations::CENTER);

  riemann = new riemannSolver(geom);

  SNES snes;
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

//  SNESSetFunction(snes, residual->globalVec, computeResidual, this);

}

timeStepper::~timeStepper()
{
  SNESDestroy(&snes);
  
  delete elem, elemOld;
  delete prim, primOld;
  delete cons, consOld;
  delete sourcesOld;
  delete residual;

  switch (params::dim)
  {
    case 1:
      delete fluxesX1;
      break;

    case 2:
      delete fluxesX1, fluxesX2;
      break;

    case 3:
      delete fluxesX1, fluxesX2, fluxesX3;
      break;
  }
}
//void timeStepper::timeStep()
//{
//  elemOld->set(primOld, geom, locations::CENTER);
//
//  elemOld->computeFluxes(geom, 0, consOld->vars);
//  elemOld->computeSources(geom, sourcesOld->vars);
////  
////  reconstruct(primOld->vars, 
////              primLeft->vars, primRight->vars,
////              directions::X1
////             );
////
//  riemann.solve(primRight->vars, primLeft->vars,
//                geom, directions::X1,
//                fluxesX1->vars
//               );
//
//}

//PetscErrorCode computeResidual(SNES snes,
//                               Vec primVec,
//                               Vec residualVec,
//                               void *ptr
//                              )
//{
//  timeStepper *ts = (class timeStepper*)ptr;
//
//}

