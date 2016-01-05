#include "reconstruction.hpp"

void reconstruction::reconstructWENO5(const grid &prim,
                                      const int dir,
                                      grid &primLeft,
                                      grid &primRight
                                     )
{
  //WENO5 algorithm, copied from SpEC (up to some left/right conventions, and
  //the use of AF...)
  const double eps2 = 1.0e-17;
  double filter1D[]  = {1., 0, 0, 0, 0, 
			0, 1., 0, 0, 0,
			0, 0, 1., 0, 0,
			0, 0, 0, 1., 0,
			0, 0, 0, 0, 1.};
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
  
  for(int var=0; var<prim.numVars; var++)
  {
  	array dvar = convolve(prim.vars[var], filter);
  	//Get stencil
  	array y0 = dvar(span,span,span,4);
  	array y1 = dvar(span,span,span,3);
  	array y2 = dvar(span,span,span,2);
  	array y3 = dvar(span,span,span,1);
  	array y4 = dvar(span,span,span,0);

  	//Compute smoothness operators
  	array beta1 = (( 4.0/3.0)*y0*y0 - (19.0/3.0)*y0*y1 +
	                 (25.0/3.0)*y1*y1 + (11.0/3.0)*y0*y2 -
        		       (31.0/3.0)*y1*y2 + (10.0/3.0)*y2*y2
                  )
                 +
                  eps2*(1.0 + af::abs(y0) + af::abs(y1) + af::abs(y2));

  	array beta2 = (( 4.0/3.0)*y1*y1 - (13.0/3.0)*y1*y2 +
	                 (13.0/3.0)*y2*y2 + ( 5.0/3.0)*y1*y3 -
        		       (13.0/3.0)*y2*y3 + ( 4.0/3.0)*y3*y3
                  ) 
                 +
              	  eps2*(1.0 + af::abs(y1) + af::abs(y2) + af::abs(y3));

	array beta3 = ((10.0/3.0)*y2*y2 - (31.0/3.0)*y2*y3 +
	                 (25.0/3.0)*y3*y3 + (11.0/3.0)*y2*y4 -
        		       (19.0/3.0)*y3*y4 + ( 4.0/3.0)*y4*y4
                  ) 
                 + 
              	  eps2*(1.0 + af::abs(y2) + af::abs(y3) + af::abs(y4));

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
  
