#include "reconstruction.hpp"

array reconstruction::minmod(array &x, array &y, array &z)
{
  array minOfAll = af::min(af::min(af::abs(x), af::abs(y)), 
		                       af::abs(z)
               			      );
  
  //Stupid convention in ArrayFire: sign(x)=1 for x<0 and sign(x)=0 for x>0
  array signx = 1.-2.*sign(x);
  array signy = 1.-2.*sign(y);
  array signz = 1.-2.*sign(z);
  
  return 0.25 * af::abs(signx + signy ) * (signx + signz ) * minOfAll;
}

array reconstruction::slopeMM(const int dir,const double dX, const array& in)
{
  double filter1D[]  = {1,-1, 0, /* Forward difference */
                			  0, 1,-1  /* Backward difference */
                       };
  array filter;
  switch (dir)
  {
    case directions::X1:
    	filter =  array(3, 1, 1, 2, filter1D)/dX;
    	break;

    case directions::X2:
    	filter =  array(1, 3, 1, 2, filter1D)/dX;
    	break;

    case directions::X3:
    	filter =  array(1, 1, 3, 2, filter1D)/dX;
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

void reconstruction::reconstructMM(const grid &prim,
                	                 const int dir,
                          		     grid &primLeft,
          		                     grid &primRight
                          		    )
{
  for(int var=0; var<prim.numVars; var++)
  {
  	//Note: we set dX=1., because the 1/dX in slope
  	//exactly cancels the dX in the computation of the
  	//value on cell faces...
  	array slope = slopeMM(dir,1.,prim.vars[var]);
  	primLeft.vars[var]  = prim.vars[var] - 0.5*slope;
  	primRight.vars[var] = prim.vars[var] + 0.5*slope;
  }
}

