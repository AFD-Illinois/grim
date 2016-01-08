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

  MinSpeedLeft = af::constant(1,
			  primLeft->vars[0].dims(directions::X1),
			  primLeft->vars[0].dims(directions::X2),
			  primLeft->vars[0].dims(directions::X3),
			  f64
			  );
  MaxSpeedLeft = af::constant(1,
			  primLeft->vars[0].dims(directions::X1),
			  primLeft->vars[0].dims(directions::X2),
			  primLeft->vars[0].dims(directions::X3),
			  f64
			  );
  MinSpeedRight = af::constant(1,
			  primLeft->vars[0].dims(directions::X1),
			  primLeft->vars[0].dims(directions::X2),
			  primLeft->vars[0].dims(directions::X3),
			  f64
			  );
  MaxSpeedRight = af::constant(1,
			  primLeft->vars[0].dims(directions::X1),
			  primLeft->vars[0].dims(directions::X2),
			  primLeft->vars[0].dims(directions::X3),
			  f64
			  );
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
  int locationLeft,locationRight;
  int shiftX1, shiftX2, shiftX3;
  int fluxDirection;
  switch (dir)
  {
    case directions::X1:
      locationLeft = locations::LEFT;
      locationRight = locations::RIGHT;
      fluxDirection = 1;
      shiftX1  = 1;
      shiftX2  = 0;
      shiftX3  = 0;
      break;

    case directions::X2:
      locationLeft = locations::BOTTOM;
      locationRight = locations::TOP;
      fluxDirection = 2;
      shiftX1  = 0;
      shiftX2  = 1;
      shiftX3  = 0;
      break;

    case directions::X3:
      locationLeft = locations::BACK;
      locationRight = locations::FRONT;
      fluxDirection = 3;
      shiftX1  = 0;
      shiftX2  = 0;
      shiftX3  = 1;
      break;
  }

  //Reconstruction gives, at a point of index i:
  // primLeft : right-biased stencil reconstructs on face i-/1.2
  // primRight: left-biased stencil reconstructs on face i+1/2
  reconstruction::reconstruct(prim, dir, *primLeft, *primRight);

  //elemLeft is the left-biased stencil on the right-face i+1/2...
  elemLeft->set(*primRight, geom, locationRight);
  //elemRight is the right-biased stencil on the left-face i-1/2...
  elemRight->set(*primLeft, geom, locationLeft);

  elemLeft->computeMinMaxCharSpeeds(geom,dir,MinSpeedLeft,MaxSpeedLeft);
  elemRight->computeMinMaxCharSpeeds(geom,dir,MinSpeedRight,MaxSpeedRight);
  
  elemLeft->computeFluxes(geom, fluxDirection, *fluxLeft);
  elemLeft->computeFluxes(geom, 0,             *consLeft);
  elemRight->computeFluxes(geom, fluxDirection, *fluxRight);
  elemRight->computeFluxes(geom, 0,             *consRight);

  for (int var=0; var<vars::dof; var++)
  {
    //The fluxes are requested on the left-face i-1/2.
    //Hence, we need to shift elemLeft,fluxLeft,consLeft by a single point to the right
    // (elemLeft[i] refers to values at i+1/2, we want values at i-1/2)
    MinSpeedLeft = af::shift(MinSpeedLeft,shiftX1, shiftX2, shiftX3);
    MaxSpeedLeft = af::shift(MaxSpeedLeft,shiftX1, shiftX2, shiftX3);
    MinSpeedLeft = min(MinSpeedLeft,MinSpeedRight);
    MaxSpeedLeft = max(MaxSpeedLeft,MaxSpeedRight);

    //HLL formula
    flux.vars[var] = (  MaxSpeedLeft*af::shift(fluxLeft->vars[var], shiftX1, shiftX2, shiftX3) 
			- MinSpeedLeft*fluxRight->vars[var]
			+ MinSpeedLeft*MaxSpeedLeft 
			* (  consRight->vars[var]
			     - af::shift(consLeft->vars[var], shiftX1, shiftX2, shiftX3)
			     )
			)/(MaxSpeedLeft-MinSpeedLeft);
  }
            
}
