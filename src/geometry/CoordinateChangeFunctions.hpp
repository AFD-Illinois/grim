
double GammieRadius(const double X1)
{
  return exp(X1);
}

double GammieTheta(const double X2)
{
  return M_PI*X2+ 0.5*(1 - params::hSlope)*sin(2.*M_PI*X2);
}

/**************************************
// The following code is largely copied
// from Sasha's version of HARM
***************************************/

//smooth step function:
// Ftr = 0 if x < 0, Ftr = 1 if x > 1 and smoothly interps. in btw.
double Ftr( const double x )
{
  double res;
  if( x <= 0 ) {
    res = 0;
  }
  else if( x >= 1 ) {
    res = 1;
  }
  else {
    res = (64 + cos(5*M_PI*x) + 70*sin((M_PI*(-1 + 2*x))/2.) + 5*sin((3*M_PI*(-1 + 2*x))/2.))/128.;
  }
  return( res );
}

double Ftrgenlin( const double x, const double xa, const double xb, const double ya, const double yb )
{
  double res = (x*ya)/xa + (-((x*ya)/xa) + ((x - xb)*(1 - yb))/(1 - xb) + yb)*Ftr((x - xa)/(-xa + xb));
  return res;
}

//goes from ya to yb as x goes from xa to xb
double Ftrgen(const double x, const double xa, const double xb, const double ya, const double yb )
{
  double res = ya + (yb-ya)*Ftr( (x-xa)/(xb-xa) );
  return res;
}

double Fangle(const double x )
{
  double res;
  if( x <= -1 ) {
    res = 0;
  }
  else if( x >= 1 ) {
    res = x;
  }
  else {
    res = (1 + x + (-140*sin((M_PI*(1 + x))/2.) + (10*sin((3*M_PI*(1 + x))/2.))/3. + (2*sin((5*M_PI*(1 + x))/2.))/5.)/(64.*M_PI))/2.;
  }
  return res;
}

double limlin( const double x, const double x0, const double dx, const double y0 )
{
  return( y0 - dx * Fangle(-(x-x0)/dx) );
}

double minlin(const double x,const double x0,const double dx,const double y0 )
{
  return( y0 + dx * Fangle((x-x0)/dx) );
}

double mins(const double f1,const double f2,const double df )
{
  return limlin(f1, f2, df, f2);
}

double maxs(const double f1,const double f2,const double df )
{
  return (-mins(-f1, -f2, df) );
}

//=mins if dir < 0
//=maxs if dir >= 0
double minmaxs(const double f1,const double f2,const double df,const double dir )
{
  if( dir>=0 ) {
    return maxs(f1, f2, df);
  }
  return mins(f1, f2, df);
}


double sinth0( double *X0, double *X)
{
  double Xc0[3];
  for(int j=0;j<3;j++)
    Xc0[j]=X[j];
  Xc0[directions::X2] = X0[directions::X2];
  return GammieRadius(X0[directions::X1]) * sin(GammieTheta(X0[directions::X2])) / GammieRadius(Xc0[directions::X1]);
}

double sinth1in( double *X0, double *X)
{
  double X0c[3];
  for(int j=0;j<3;j++)
    X0c[j]=X0[j];
  X0c[directions::X2] = X[directions::X2];
  return GammieRadius(X0[directions::X1]) * sin(GammieTheta(X0c[directions::X2])) / GammieRadius(X[directions::X1]);
}


double th2in( double *X0, double *X)
{
  double Xc0[3];
  double Xcmid[3];
  double res;
  double th0;
  for(int j=0;j<3;j++)
    Xc0[j]=X[j];
  Xc0[directions::X2] = X0[directions::X2];
  for(int j=0;j<3;j++)
    Xcmid[j]=X[j];
  Xcmid[directions::X2] = 0.5;
  
  th0 = asin( sinth0(X0, X) );
  res = (GammieTheta(X[directions::X2]) - GammieTheta(Xc0[directions::X2]))/(GammieTheta(Xcmid[directions::X2]) - GammieTheta(Xc0[directions::X2])) * 
    (GammieTheta(Xcmid[directions::X2])-th0) + th0;
  return( res );
}

double func1( double *X0, double *X)
{
  return sin(GammieTheta(X[directions::X2]));
}

double func2( double *X0, double *X)
{
  double Xca[3];
  double func2;
  double sth1in, sth2in, sth1inaxis, sth2inaxis;
 
  for(int j=0;j<3;j++)
    Xca[j]=X[j];
  Xca[directions::X2] = 0.;
  
  sth1in = sinth1in( X0, X);
  sth2in = sin( th2in(X0, X));
  
  sth1inaxis = sinth1in( X0, Xca);
  sth2inaxis = sin( th2in(X0, Xca) );
  
  func2 = minmaxs( sth1in, sth2in, fabs(sth2inaxis-sth1inaxis)+1.e-20, X[directions::X1] - X0[directions::X1] );
  
  return( func2 );
}
