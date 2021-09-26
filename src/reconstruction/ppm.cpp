#include "reconstruction.hpp"

array reconstruction::slopePPM(const int dir,const double dX, const array& in,
			       int &numReads, int &numWrites
			       )
{
 // Construct stencil
  int x0Shift, y0Shift, z0Shift;
  int x1Shift, y1Shift, z1Shift;
  int x3Shift, y3Shift, z3Shift;
  int x4Shift, y4Shift, z4Shift;
  switch (dir)
  {
    case directions::X1:
      x0Shift = 2; y0Shift = 0;  z0Shift = 0;
      x1Shift = 1; y1Shift = 0;  z1Shift = 0;
      x3Shift = -1; y3Shift = 0;  z3Shift = 0;
      x4Shift = -2; y4Shift = 0;  z4Shift = 0;
      break;

    case directions::X2:
      x0Shift =  0; y0Shift = 2; z0Shift = 0;
      x1Shift =  0; y1Shift = 1; z1Shift = 0;
      x3Shift =  0; y3Shift = -1; z3Shift = 0;
      x4Shift =  0; y4Shift = -2; z4Shift = 0;
      break;

    case directions::X3:
      x0Shift =  0; y0Shift =  0; z0Shift = 2;
      x1Shift =  0; y1Shift =  0; z1Shift = 1;
      x3Shift =  0; y3Shift =  0; z3Shift = -1;
      x4Shift =  0; y4Shift =  0; z4Shift = -2;
      break;
  }
  
  array y0 = af::shift(in, x0Shift, y0Shift, z0Shift);
  array y1 = af::shift(in, x1Shift, y1Shift, z1Shift);
  array y2 = in;
  array y3 = af::shift(in, x3Shift, y3Shift, z3Shift);
  array y4 = af::shift(in, x4Shift, y4Shift, z4Shift);
    
  // Approximants for slopes
  array d0 = 2.*(y1-y0);
  array d1 = 2.*(y2-y1);
  array d2 = 2.*(y3-y2);
  array d3 = 2.*(y4-y3);   
  array D1 = 0.5*(y2-y0);
  array D2 = 0.5*(y3-y1);
  array D3 = 0.5*(y4-y2);
  
  array condZeroSlope1 = (d1*d0<=0.);
  array sign1 = (D1>0.)*2.-1.;
  array DQ1 = (1.-condZeroSlope1)*sign1*af::min(af::abs(D1),af::min(af::abs(d0),af::abs(d1)));
  array condZeroSlope2 = (d2*d1<=0.);
  array sign2 = (D2>0.)*2.-1.;
  array DQ2 = (1.-condZeroSlope2)*sign2*af::min(af::abs(D2),af::min(af::abs(d1),af::abs(d2)));
  array condZeroSlope3 = (d3*d2<=0.);
  array sign3 = (D3>0.)*2.-1.;
  array DQ3 = (1.-condZeroSlope3)*sign3*af::min(af::abs(D3),af::min(af::abs(d2),af::abs(d3)));
  
  // Base high-order PPM reconstruction
  array leftV = 0.5*(y2+y1)-1./6.*(DQ2-DQ1);
  array rightV = 0.5*(y3+y2)-1./6.*(DQ3-DQ2);
  
  // Corrections
  array corr1 = ((rightV-y2)*(y2-leftV)<=0.);
  array qd = rightV-leftV;
  array qe = 6.*(y2-0.5*(rightV+leftV));
  array corr2 = (qd*(qd-qe)<0.);
  array corr3 = (qd*(qd+qe)<0.);
  leftV = leftV*(1.-corr1)+corr1*y2;
  rightV = rightV*(1.-corr1)+corr1*y2;
  
  leftV = leftV*(1.-corr2)+corr2*(3.*y2-2.*rightV);
  rightV = rightV*corr2+(1.-corr2)*rightV*(1.-corr3)+(1.-corr2)*corr3*(3.*y2-2.*leftV);
  return (rightV-leftV)/dX;
}

// PPM algorithm, adapted from HARM code (by X. Guan)
// ref. Colella && Woodward's PPM paper
void reconstruction::reconstructPPM(const grid &prim,
				    const int dir,
				    grid &primLeft,
				    grid &primRight,
				    int &numReads,
                                    int &numWrites
				    )
{
  // Construct stencil
  int x0Shift, y0Shift, z0Shift;
  int x1Shift, y1Shift, z1Shift;
  int x3Shift, y3Shift, z3Shift;
  int x4Shift, y4Shift, z4Shift;
  switch (dir)
  {
    case directions::X1:
      x0Shift = 2; y0Shift = 0;  z0Shift = 0;
      x1Shift = 1; y1Shift = 0;  z1Shift = 0;
      x3Shift = -1; y3Shift = 0;  z3Shift = 0;
      x4Shift = -2; y4Shift = 0;  z4Shift = 0;
    	break;

    case directions::X2:
      x0Shift =  0; y0Shift = 2; z0Shift = 0;
      x1Shift =  0; y1Shift = 1; z1Shift = 0;
      x3Shift =  0; y3Shift = -1; z3Shift = 0;
      x4Shift =  0; y4Shift = -2; z4Shift = 0;
    	break;

    case directions::X3:
      x0Shift =  0; y0Shift =  0; z0Shift = 2;
      x1Shift =  0; y1Shift =  0; z1Shift = 1;
      x3Shift =  0; y3Shift =  0; z3Shift = -1;
      x4Shift =  0; y4Shift =  0; z4Shift = -2;
    	break;
  }
  
  for(int var=0; var<prim.numVars; var++)
  {
    array y0 = af::shift(prim.vars[var], x0Shift, y0Shift, z0Shift);
    array y1 = af::shift(prim.vars[var], x1Shift, y1Shift, z1Shift);
    array y2 = prim.vars[var];
    array y3 = af::shift(prim.vars[var], x3Shift, y3Shift, z3Shift);
    array y4 = af::shift(prim.vars[var], x4Shift, y4Shift, z4Shift);
    
    // Approximants for slopes
    array d0 = 2.*(y1-y0);
    array d1 = 2.*(y2-y1);
    array d2 = 2.*(y3-y2);
    array d3 = 2.*(y4-y3);   
    array D1 = 0.5*(y2-y0);
    array D2 = 0.5*(y3-y1);
    array D3 = 0.5*(y4-y2);
    
    array condZeroSlope1 = (d1*d0<=0.);
    array sign1 = (D1>0.)*2.-1.;
    array DQ1 = (1.-condZeroSlope1)*sign1*af::min(af::abs(D1),af::min(af::abs(d0),af::abs(d1)));
    array condZeroSlope2 = (d2*d1<=0.);
    array sign2 = (D2>0.)*2.-1.;
    array DQ2 = (1.-condZeroSlope2)*sign2*af::min(af::abs(D2),af::min(af::abs(d1),af::abs(d2)));
    array condZeroSlope3 = (d3*d2<=0.);
    array sign3 = (D3>0.)*2.-1.;
    array DQ3 = (1.-condZeroSlope3)*sign3*af::min(af::abs(D3),af::min(af::abs(d2),af::abs(d3)));

    // Base high-order PPM reconstruction
    array leftV = 0.5*(y2+y1)-1./6.*(DQ2-DQ1);
    array rightV = 0.5*(y3+y2)-1./6.*(DQ3-DQ2);
    
    // Corrections
    array corr1 = ((rightV-y2)*(y2-leftV)<=0.);
    array qd = rightV-leftV;
    array qe = 6.*(y2-0.5*(rightV+leftV));
    array corr2 = (qd*(qd-qe)<0.);
    array corr3 = (qd*(qd+qe)<0.);
    leftV = leftV*(1.-corr1)+corr1*y2;
    rightV = rightV*(1.-corr1)+corr1*y2;
    
    leftV = leftV*(1.-corr2)+corr2*(3.*y2-2.*rightV);
    rightV = rightV*corr2+(1.-corr2)*rightV*(1.-corr3)+(1.-corr2)*corr3*(3.*y2-2.*leftV);

    primLeft.vars[var] = leftV;
    primRight.vars[var] = rightV;
  }
}
