#include "reconstruction.hpp"

array reconstruction::slopeWENO5(const int dir,const double dX, const array& in,
                                 int &numReads, int &numWrites
                                )
{
  //WENO5 algorithm, copied from SpEC (up to some left/right conventions, and
  //the use of AF...)
  const double eps2 = 1.0e-17;

  int x0Shift, y0Shift, z0Shift;
  int x1Shift, y1Shift, z1Shift;
  int x3Shift, y3Shift, z3Shift;
  int x4Shift, y4Shift, z4Shift;
  switch (dir)
  {
    case directions::X1:
      x0Shift = -2; y0Shift = 0;  z0Shift = 0;
      x1Shift = -1; y1Shift = 0;  z1Shift = 0;
      x3Shift =  1; y3Shift = 0;  z3Shift = 0;
      x4Shift =  2; y4Shift = 0;  z4Shift = 0;
    	break;

    case directions::X2:
      x0Shift =  0; y0Shift = -2; z0Shift = 0;
      x1Shift =  0; y1Shift = -1; z1Shift = 0;
      x3Shift =  0; y3Shift =  1; z3Shift = 0;
      x4Shift =  0; y4Shift =  2; z4Shift = 0;
    	break;

    case directions::X3:
      x0Shift =  0; y0Shift =  0; z0Shift = -2;
      x1Shift =  0; y1Shift =  0; z1Shift = -1;
      x3Shift =  0; y3Shift =  0; z3Shift =  1;
      x4Shift =  0; y4Shift =  0; z4Shift =  2;
    	break;
  }

  array y0 = af::shift(in, x0Shift, y0Shift, z0Shift);
  array y1 = af::shift(in, x1Shift, y1Shift, z1Shift);
  array y2 = in;
  array y3 = af::shift(in, x3Shift, y3Shift, z3Shift);
  array y4 = af::shift(in, x4Shift, y4Shift, z4Shift);
  /* Reads:
   * -----
   * shift -2, -1, 1, 2: 4
   *
   * Writes:
   * ------
   * y0, y1, y3, y4 : 4 */

  //temporary formula : smooth 5pt stencil
  array ans = (-y4+8.*y3-8.*y1+y0)/12./dX;
  ans.eval();
  /* Reads:
   * -----
   *  y4, y3, y1, y0 : 4
   *
   * Writes:
   * ------
   * ans : 1 */
  numReads  = 8;
  numWrites = 5;
  
  return ans;

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
  
  //TODO: Compute weights
  array w1r = 1.0/(16.0*beta1*beta1);
  array w2r = 5.0/( 8.0*beta2*beta2);
  array w3r = 5.0/(16.0*beta3*beta3);
  array w1l = 5.0/(16.0*beta1*beta1);
  array w2l = 5.0/( 8.0*beta2*beta2);
  array w3l = 1.0/(16.0*beta3*beta3);
  array denl = w1l + w2l + w3l;
  array denr = w1r + w2r + w3r;
  
  //TODO: Substencil Interpolations
  array u1r =  0.375*y0 - 1.25*y1 + 1.875*y2;
  array u2r = -0.125*y1 + 0.75*y2 + 0.375*y3;
  array u3r =  0.375*y2 + 0.75*y3 - 0.125*y4;
  array u1l = -0.125*y0 + 0.75*y1 + 0.375*y2;
  array u2l =  0.375*y1 + 0.75*y2 - 0.125*y3;
  array u3l =  1.875*y2 - 1.25*y3 + 0.375*y4;
  
  //Reconstruction
  
}

void reconstruction::reconstructWENO5(const grid &prim,
                                      const int dir,
                                      grid &primLeft,
                                      grid &primRight,
                                      int &numReads,
                                      int &numWrites
                                     )
{
  //WENO5 algorithm, copied from SpEC (up to some left/right conventions, and
  //the use of AF...)
  const double eps2 = 1.0e-17;

  int x0Shift, y0Shift, z0Shift;
  int x1Shift, y1Shift, z1Shift;
  int x3Shift, y3Shift, z3Shift;
  int x4Shift, y4Shift, z4Shift;
  switch (dir)
  {
    case directions::X1:
      x0Shift = -2; y0Shift = 0;  z0Shift = 0;
      x1Shift = -1; y1Shift = 0;  z1Shift = 0;
      x3Shift =  1; y3Shift = 0;  z3Shift = 0;
      x4Shift =  2; y4Shift = 0;  z4Shift = 0;
    	break;

    case directions::X2:
      x0Shift =  0; y0Shift = -2; z0Shift = 0;
      x1Shift =  0; y1Shift = -1; z1Shift = 0;
      x3Shift =  0; y3Shift =  1; z3Shift = 0;
      x4Shift =  0; y4Shift =  2; z4Shift = 0;
    	break;

    case directions::X3:
      x0Shift =  0; y0Shift =  0; z0Shift = -2;
      x1Shift =  0; y1Shift =  0; z1Shift = -1;
      x3Shift =  0; y3Shift =  0; z3Shift =  1;
      x4Shift =  0; y4Shift =  0; z4Shift =  2;
    	break;
  }
  
  for(int var=0; var<prim.numVars; var++)
  {
    array y0 = af::shift(prim.vars[var], x0Shift, y0Shift, z0Shift);
    array y1 = af::shift(prim.vars[var], x1Shift, y1Shift, z1Shift);
    array y2 = prim.vars[var];
    array y3 = af::shift(prim.vars[var], x3Shift, y3Shift, z3Shift);
    array y4 = af::shift(prim.vars[var], x4Shift, y4Shift, z4Shift);
    /* Reads:
     * -----
     * shift -2, -1, 1, 2: 4*prim.numVars
     *
     * Writes:
     * ------
     * y0, y1, y3, y4 : 4*prim.numvars */

  	//Compute smoothness operators
  	array beta1 = (( 4.0/3.0)*y0*y0 - (19.0/3.0)*y0*y1 +
	                 (25.0/3.0)*y1*y1 + (11.0/3.0)*y0*y2 -
        		       (31.0/3.0)*y1*y2 + (10.0/3.0)*y2*y2
                  )
                 +
                  eps2*(1.0 + af::abs(y0) + af::abs(y1) + af::abs(y2));
    beta1.eval();
    /* Reads:
     * -----
     * y0, y1, y2: 3*prim.numVars
     *
     * Writes:
     * ------
     * beta1: prim.numvars */

  	array beta2 = (( 4.0/3.0)*y1*y1 - (13.0/3.0)*y1*y2 +
	                 (13.0/3.0)*y2*y2 + ( 5.0/3.0)*y1*y3 -
        		       (13.0/3.0)*y2*y3 + ( 4.0/3.0)*y3*y3
                  ) 
                 +
              	  eps2*(1.0 + af::abs(y1) + af::abs(y2) + af::abs(y3));
    beta2.eval();
    /* Reads:
     * -----
     * y1, y2, y3: 3*prim.numVars
     *
     * Writes:
     * ------
     * beta2: prim.numvars */

	  array beta3 = ((10.0/3.0)*y2*y2 - (31.0/3.0)*y2*y3 +
	                 (25.0/3.0)*y3*y3 + (11.0/3.0)*y2*y4 -
        		       (19.0/3.0)*y3*y4 + ( 4.0/3.0)*y4*y4
                  ) 
                 + 
              	  eps2*(1.0 + af::abs(y2) + af::abs(y3) + af::abs(y4));
    beta3.eval();
    /* Reads:
     * -----
     * y2, y3, y4: 3*prim.numVars
     *
     * Writes:
     * ------
     * beta3: prim.numvars */

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

    primLeft.vars[var].eval();
    /* Reads:
     * -----
     * beta1, beta2, beta3, y0, y1, y2, y3, y4: 7*prim.numVars
     *
     * Writes:
     * ------
     * primLeft: prim.numvars */

    primRight.vars[var].eval();
    /* Reads:
     * -----
     * beta1, beta2, beta3, y0, y1, y2, y3, y4: 7*prim.numVars
     *
     * Writes:
     * ------
     * primRight: prim.numvars */
  }
  /* Total reads : 27*prim.numVars
   * Total writes: 9*prim.numVars */ 

  numReads  = 27*prim.numVars;
  numWrites = 9*prim.numVars;
}
  
