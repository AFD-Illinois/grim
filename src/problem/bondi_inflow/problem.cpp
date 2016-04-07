#include "../problem.hpp"

namespace vars
{
  int dof = 8;
  int Q = 8;
  int DP = 9;
};


namespace params
{
  int N1 = 20;
  int N2 = 20;
  int N3 = 1;
  int dim = 2;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double InitialDt = .2/N1;
  double Time = 0.;
  double finalTime = 0.5;
  int metric = metrics::MODIFIED_KERR_SCHILD;
  double hSlope = 1.0;
  double blackHoleSpin = 0.0;

  double sonicRadius = 8.;
  double mDot = 2.;
  double bMag = 1.;

  double X1Start = log(3.1), X1End = log(35.);
  double X2Start = 0.+1.e-8, X2End = 1.-1.e-8;
  double X3Start = 0., X3End = M_PI/10.;

  int boundaryLeft   = boundaries::DIRICHLET;
  int boundaryRight  = boundaries::DIRICHLET;

  int boundaryTop    = boundaries::DIRICHLET;
  int boundaryBottom = boundaries::DIRICHLET;

  int boundaryFront  = boundaries::DIRICHLET;
  int boundaryBack   = boundaries::DIRICHLET;

  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;

  int conduction = 0;
  int viscosity  = 0;
  int highOrderTermsConduction = 0.;
  int highOrderTermsViscosity = 0.;
  double adiabaticIndex = 4./3;

  double slopeLimTheta = 2;
  int reconstruction = reconstructionOptions::WENO5;
  int riemannSolver  = riemannSolvers::HLL;

  int maxNonLinearIter = 10;
  int maxLineSearchIters = 10;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-16;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;
  
};

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
}


double FuncT(const double T, const double R, const double C1, const double C2)
{
  double nPoly = 1./(params::adiabaticIndex-1.);
  return pow(1.+(1.+nPoly)*T,2.)*(1.-2./R+pow(C1/R/R/pow(T,nPoly),2.))-C2;
}

double SolveT(const double R, const double C1, const double C2)
{
  double nPoly = 1./(params::adiabaticIndex-1.);
  double rtol = 1.e-12;
  double ftol = 1.e-14;
  double Tmin = 0.6*(sqrt(C2)-1.)/(nPoly+1.);
  double Tmax = pow(C1*sqrt(2./R/R/R),1./nPoly);
  double f0,f1,fh;
  double T0,T1,Th;
  T0=0.6*Tmin;
  f0=FuncT(T0,R,C1,C2);
  T1=Tmax;
  f1=FuncT(T1,R,C1,C2);
  if(f0*f1>0.)
    {
      printf("Failed solving for T at R = %f; C1 = %f; C2 = %f \n",R,C1,C2);
      PetscPrintf(PETSC_COMM_WORLD, "Failed determination of T \n");
      exit(1);
    }
  Th = (f1*T0-f0*T1)/(f1-f0);
  fh = FuncT(Th,R,C1,C2);
  double EpsT = rtol*(Tmin+Tmax);
  while(fabs(Th-T0)>EpsT && fabs(Th-T1)>EpsT && fabs(fh)>ftol)
    {
      if(fh*f0<0.)
	{
	  T0=Th;
	  f0=fh;
	}
      else
	{
	  T1=Th;
	  f1=fh;
	}
      Th = (f1*T0-f0*T1)/(f1-f0);
      fh = FuncT(Th,R,C1,C2);
    }
  return Th;
}

void timeStepper::initialConditions(int &numReads,int &numWrites)
{
  const double nPoly = 1./(params::adiabaticIndex-1.);
  const double Rc = params::sonicRadius; 
  const double uc = sqrt(0.5/Rc);
  const double vc2 = nPoly*uc*uc/(1.-3.*uc*uc);
  const double Tc = vc2/(1.-vc2)/(nPoly+1.);
  const double C1 = pow(Tc,nPoly)*uc*Rc*Rc;
  const double C2 = pow(1.+(nPoly+1.)*Tc,2.)*(1.-2./Rc+uc*uc);
  const double Kp = pow(params::mDot/4./M_PI/C1,-1./nPoly);

  array xCoords[3];
  for(int d=0;d<3;d++)
    xCoords[d]=XCoords->vars[d];
  geomCenter->XCoordsToxCoords(XCoords->vars,xCoords);

  const array& R = xCoords[0];

  //Ok, let's ignore ArrayFire here - it's only the initial
  //conditions of a test...
  const int N1g = params::N1+2*params::numGhost;
  const int N2g = params::dim>1 ? params::N2+2*params::numGhost : 1;
  const int N3g = params::dim>2 ? params::N3+2*params::numGhost : 1;

  printf("Using grid of %i x %i x %i\n",N1g,N2g,N3g);
 
  double* host_R = R.host<double>();
  double T[N1g];
  double Rho0[N1g];
  double ur[N1g];
  double ut[N1g];

  for(int i=0;i<N1g;i++)
    {
      const double Rl = host_R[i];
      double Tl,rhol,url,utl;
      if(Rl>2.+1.e-8)
	{
	  Tl = SolveT(Rl,C1,C2);
	  rhol = pow(Tl/Kp,nPoly);
	  url = -C1/Rl/Rl/pow(Tl,nPoly);
	  utl = (2.*url/Rl+sqrt(url*url+1.-2./Rl))/(1.-2./Rl);
	}
      else
	{
	  Tl = 1.; rhol=params::rhoFloorInFluidElement;
	  url=0.; utl=1.;
	}
      T[i]=Tl; Rho0[i]=rhol; ur[i]=url; ut[i]=utl; 
      //printf("r=%e; rho=%e; T=%e;\n",Rl,rhol,Tl);
    }

  //Set initial conditions
  primOld->vars[vars::RHO] = 1.;
  primOld->vars[vars::U] = 0.;
  primOld->vars[vars::U1] = 0.;
  primOld->vars[vars::B1] = 0.;
  for(int i=0;i<N1g;i++)
    {
      const double Rl = host_R[i];
      for(int j=0;j<N2g;j++)
	for(int k=0;k<N3g;k++)
	  {
	    primOld->vars[vars::RHO](i,j,k) = Rho0[i];
	    primOld->vars[vars::U](i,j,k) = nPoly*Rho0[i]*T[i];
	    primOld->vars[vars::U1](i,j,k) = (ur[i] + ut[i]*2./(Rl+2.))/Rl;
	    primOld->vars[vars::B1](i,j,k)  = params::bMag/Rl/Rl/Rl; 
	  }
    }
  primOld->vars[vars::U2]  = 0.; 
  primOld->vars[vars::U3]  = 0.; 
  primOld->vars[vars::B2]  = 0.;
  primOld->vars[vars::B3]  = 0.;
  
  for (int var=0; var<vars::dof; var++)
  {
    primOld->vars[var].eval();
  }

  for (int var=0; var<vars::dof; var++) 
    {
      primIC->vars[var]=1.*primOld->vars[var];
      primIC->vars[var].eval();
    }

  af::sync();
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{
  printf("halfStepDiagnostics\n");
}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
  printf("fullStepDiagnostics: Time = %e\n",params::Time);
  //af_print(primOld->varsOld->vars[vars::RHO],8.);
}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{
  for (int var=0; var<vars::dof; var++) 
    {
      for(int i=0;i<params::numGhost;i++)
	primOld->vars[var](i,span,span)=primIC->vars[var](i,span,span);
      for(int i=params::N1+params::numGhost;i<params::N1+2*params::numGhost;i++)
	primOld->vars[var](i,span,span)=primIC->vars[var](i,span,span);
      if(params::dim==1)
	continue;
      for(int i=0;i<params::numGhost;i++)
	primOld->vars[var](span,i,span)=primIC->vars[var](span,i,span);
      for(int i=params::N2+params::numGhost;i<params::N2+2*params::numGhost;i++)
	primOld->vars[var](span,i,span)=primIC->vars[var](span,i,span);
      if(params::dim==2)
	continue;
      for(int i=0;i<params::numGhost;i++)
	primOld->vars[var](span,span,i)=primIC->vars[var](span,span,i);
      for(int i=params::N3+params::numGhost;i<params::N3+2*params::numGhost;i++)
	primOld->vars[var](span,span,i)=primIC->vars[var](span,span,i);
    }
  for (int var=0; var<vars::dof; var++)
    primOld->vars[var].eval();
  af::sync();
};
