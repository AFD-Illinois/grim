#include "reconstruction.hpp"

array reconstruction::minmod(array &x, array &y, array &z,
                             int &numReads, int &numWrites
                            )
{
  array minOfAll = af::min(af::min(af::abs(x), af::abs(y)), 
		                       af::abs(z)
               			      );
  
  //Stupid convention in ArrayFire: sign(x)=1 for x<0 and sign(x)=0 for x>0
  array signx = 1.-2.*sign(x);
  array signy = 1.-2.*sign(y);
  array signz = 1.-2.*sign(z);
  
  array result = 0.25 * af::abs(signx + signy ) * (signx + signz ) * minOfAll;

  return result;
}

array reconstruction::slopeMM(const int dir,const double dX, const array& in,
                              int &numReads,
                              int &numWrites
                             )
{
  int X1_plus,  X2_plus,  X3_plus;
  int X1_minus, X2_minus, X3_minus;
  switch (dir)
  {
    case directions::X1:
	
	X1_plus  = -1; X2_plus  = 0; X3_plus  = 0;
	X1_minus =  1; X2_minus = 0; X3_minus = 0;

    	break;

    case directions::X2:
	
	X1_plus  = 0; X2_plus  = -1; X3_plus  = 0;
	X1_minus = 0; X2_minus =  1; X3_minus = 0;

    	break;

    case directions::X3:
	
	X1_plus  = 0; X2_plus  = 0; X3_plus  = -1;
	X1_minus = 0; X2_minus = 0; X3_minus =  1;

    	break;
  }
  
  array in_plus  = af::shift(in, X1_plus,  X2_plus,  X3_plus );
  array in_minus = af::shift(in, X1_minus, X2_minus, X3_minus);

  array forwardDiff  = (in_plus -  in)/dX;
  array backwardDiff = (in - in_minus)/dX;
  array centralDiff  = backwardDiff + forwardDiff;
  
  /* TODO: add an argument to slopeLimTheta.*/
  double slopeLimTheta = params::slopeLimTheta;
  array left   = slopeLimTheta * backwardDiff;
  array center = 0.5 * centralDiff;
  array right  = slopeLimTheta * forwardDiff;
  
  int numReadsMinMod, numWritesMinMod;
  array result = minmod(left, center, right,
                        numReadsMinMod, numWritesMinMod
                       );

  return result;
}

void reconstruction::reconstructMM(const grid &prim,
                	                 const int dir,
                          		     grid &primLeft,
          		                     grid &primRight,
                                   int &numReads,
                                   int &numWrites
                          		    )
{
  int numReadsMinMod, numWritesMinMod;
  for(int var=0; var<prim.numVars; var++)
  {
  	//Note: we set dX=1., because the 1/dX in slope
  	//exactly cancels the dX in the computation of the
  	//value on cell faces...
  	array slope = slopeMM(dir,1.,prim.vars[var],
                          numReadsMinMod, numWritesMinMod
                         );
  	primLeft.vars[var]  = prim.vars[var] - 0.5*slope;
  	primRight.vars[var] = prim.vars[var] + 0.5*slope;
	
	//af::eval(primLeft.vars[var], primRight.vars[var]);
  }
}

