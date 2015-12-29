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
      shiftX1  = -1;
      shiftX2  = 0;
      shiftX3  = 0;
      break;

    case directions::X2:
      location = locations::BOTTOM;
      fluxDirection = 2;
      shiftX1  = 0;
      shiftX2  = -1;
      shiftX3  = 0;
      break;

    case directions::X3:
      location = locations::BACK;
      fluxDirection = 3;
      shiftX1  = 0;
      shiftX2  = 0;
      shiftX3  = -1;
      break;
  }

  reconstruct(prim, dir, *primLeft, *primRight);

  elemLeft->set(*primRight, geom, location);
  elemRight->set(*primLeft, geom, location);

  elemLeft->computeFluxes(geom, fluxDirection, *fluxLeft);
  elemLeft->computeFluxes(geom, 0,             *consLeft);

  elemRight->computeFluxes(geom, fluxDirection, *fluxRight);
  elemRight->computeFluxes(geom, 0,             *consRight);

  double cLaxFriedrichs = 1.;

  for (int var=0; var<vars::dof; var++)
  {
    flux.vars[var] = 0.5*(  fluxLeft->vars[var] 
                          + af::shift(fluxRight->vars[var], shiftX1, shiftX2, shiftX3)
                          - cLaxFriedrichs 
                          * (  af::shift(consRight->vars[var], shiftX1, shiftX2, shiftX3)
                             - consLeft->vars[var]
                            )
                         );
  }
            
}

array minmod(array &x, array &y, array &z)
{
  array minOfAll = af::min(af::min(af::abs(x), af::abs(y)), 
                           af::abs(z)
                          );

  //Stupid convention in ArrayFire: sign(x)=1 for x<0 and sign(x)=0 for x>0
  array signx = 1.-2.*sign(x);
  array signy = 1.-2.*sign(y);
  array signz = 1.-2.*sign(z);

  return 0.25 * abs(signx + signy ) * (signx + signz ) * minOfAll;
}

array riemannSolver::slopeMM(const int dir, const array& in)
{
  double filter1D[]  = {1,-1, 0, /* Forward difference */
                        0, 1,-1  /* Backward difference */
                       };
  array filter;
  switch (dir)
  {
    case directions::X1:
      filter =  array(3, 1, 1, 2, filter1D)/primLeft->dX1;
      break;

    case directions::X2:
      filter =  array(1, 3, 1, 2, filter1D)/primLeft->dX2;
      break;

    case directions::X3:
      filter =  array(1, 1, 3, 2, filter1D)/primLeft->dX3;
      break;
  }

  array dvar_dX = convolve(in, filter);

  array forwardDiff  = dvar_dX(span, span, span, 0);
  array backwardDiff = dvar_dX(span, span, span, 1);
  array centralDiff  = backwardDiff + forwardDiff;

  array left   = params::slopeLimTheta * backwardDiff;
  array center = 0.5 * centralDiff;
  array right  = params::slopeLimTheta * forwardDiff;

  return  minmod(left, center, right);
}

void riemannSolver::reconstructWENO5(const grid &prim,
				     const int dir,
				     grid &primLeft,
				     grid &primRight
				     )
{
  //WENO5 algorithm, copied from SpEC (up to some left/right conventions, and the use of AF...)
  const double eps2 = 1.0e-17;
  double filter1D[]  = {-1, 0, 0, 0, 0, 
			0, -1, 0, 0, 0,
			0, 0, -1, 0, 0,
			0, 0, 0, -1, 0,
			0, 0, 0, 0, -1};
  array filter;
  switch (dir)
    {
    case directions::X1:
      filter =  array(5, 1, 1, 5, filter1D);
      break;

    case directions::X2:
      filter =  array(1, 5, 1, 5, filter1D);
      break;

    case directions::X3:
      filter =  array(1, 1, 5, 5, filter1D);
      break;
    }

  for(int var=0; var<vars::dof; var++)
    {
      array dvar = convolve(prim.vars[var], filter);
      //Get stencil
      array y0 = dvar(span,span,span,0);
      array y1 = dvar(span,span,span,1);
      array y2 = dvar(span,span,span,2);
      array y3 = dvar(span,span,span,3);
      array y4 = dvar(span,span,span,4);

      //Compute smoothness operators
      array beta1 = (( 4.0/3.0)*y0*y0 - (19.0/3.0)*y0*y1 +
		     (25.0/3.0)*y1*y1 + (11.0/3.0)*y0*y2 -
		     (31.0/3.0)*y1*y2 + (10.0/3.0)*y2*y2) +
	eps2*(1.0 + abs(y0) + abs(y1) + abs(y2));
      array beta2 = (( 4.0/3.0)*y1*y1 - (13.0/3.0)*y1*y2 +
		     (13.0/3.0)*y2*y2 + ( 5.0/3.0)*y1*y3 -
		     (13.0/3.0)*y2*y3 + ( 4.0/3.0)*y3*y3) +                                                                                                    
        eps2*(1.0 + abs(y1) + abs(y2) + abs(y3));
      array beta3 = ((10.0/3.0)*y2*y2 - (31.0/3.0)*y2*y3 +
		      (25.0/3.0)*y3*y3 + (11.0/3.0)*y2*y4 -
		      (19.0/3.0)*y3*y4 + ( 4.0/3.0)*y4*y4) + 
        eps2*(1.0 + abs(y2) + abs(y3) + abs(y4));

      //Compute weights
      array w1r = 1.0/(16.0*beta1*beta1);
      array w2r = 5.0/( 8.0*beta2*beta2);
      array w3r = 5.0/(16.0*beta3*beta3);
      array w1l = 5.0/(16.0*beta1*beta1);
      array w2l = 5.0/( 8.0*beta2*beta2);
      array w3l = 1.0/(16.0*beta3*beta3);
      array denl = w1l + w2l + w3l;
      array denr = w1r + w2r + w3r;

      // Substencil Interpolations
      array u1r =  0.375*y0 - 1.25*y1 + 1.875*y2;
      array u2r = -0.125*y1 + 0.75*y2 + 0.375*y3;
      array u3r =  0.375*y2 + 0.75*y3 - 0.125*y4;
      array u1l = -0.125*y0 + 0.75*y1 + 0.375*y2;
      array u2l =  0.375*y1 + 0.75*y2 - 0.125*y3;
      array u3l =  1.875*y2 - 1.25*y3 + 0.375*y4;

      //Reconstruction
      primLeft.vars[var] = (w1l*u1l + w2l*u2l + w3l*u3l) / denl;
      primRight.vars[var] = (w1r*u1r + w2r*u2r + w3r*u3r) / denr;
    }
}


void riemannSolver::reconstruct(const grid &prim,
                                const int dir,
                                grid &primLeft,
                                grid &primRight
                               )
{
  for(int var=0; var<vars::dof; var++)
  {
    array slope = slopeMM(dir,prim.vars[var]);

    primLeft.vars[var]  = prim.vars[var] - 0.5*slope;
    primRight.vars[var] = prim.vars[var] + 0.5*slope;
  }

}
