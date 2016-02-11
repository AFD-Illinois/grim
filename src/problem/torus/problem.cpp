#include "torus.hpp"

namespace vars
{
  int dof = 8;
  int Q = 8;
  int DP = 9;
};


namespace params
{
  int N1 = 128;
  int N2 = 128;
  int N3 = 1;
  int dim = 2;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double dt = .01;
  double Time = 0.;
  double finalTime = 0.5;
  int metric = metrics::MODIFIED_KERR_SCHILD;
  double hSlope = 0.3;

  double blackHoleSpin = 0.9375;
  double InnerEdgeRadius = 6.;
  double PressureMaxRadius = 12.;
  double MinPlasmaBeta  = 15.;
  double MagneticLoops = 1;
  double Adiabat = 0.001;

  double Rin = 0.98*(1.+sqrt(1.-blackHoleSpin*blackHoleSpin));
  double Rout = 40.;
  
  double X1Start = log(Rin), X1End = log(Rout);
  double X2Start = 0.+1.e-8, X2End = 1.-1.e-8;
  double X3Start = 0., X3End = M_PI*2.;

  int boundaryLeft   = boundaries::OUTFLOW;
  int boundaryRight  = boundaries::OUTFLOW;

  int boundaryTop    = boundaries::OUTFLOW;
  int boundaryBottom = boundaries::OUTFLOW;

  int boundaryFront  = boundaries::PERIODIC;
  int boundaryBack   = boundaries::PERIODIC;

  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;

  int conduction = 0;
  int viscosity  = 0;
  int highOrderTermsConduction = 1.;
  int highOrderTermsViscosity = 1.;
  double adiabaticIndex = 4./3;

  double slopeLimTheta = 2;
  int reconstruction = reconstructionOptions::WENO5;
  int riemannSolver  = riemannSolvers::HLL;

  int maxNonLinearIter = 10;
  int maxLineSearchIters = 10;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-10;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;
  
};

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
}


/* Calculate the constant angular momentum per unit inertial mass (l = u_phi *
 * u^t) for a given black hole spin and a radius of the accretion disk.  Eqn 3.8
 * of Fishbone and Moncrief, 1976 */
double lFishboneMoncrief(double a, double r, double theta)
{
  double M = 1.;
  return sqrt(M/pow(r, 3.)) \
        *(  pow(r, 4.) + r*r*a*a - 2.*M*r*a*a \
          - a*sqrt(M*r)*(r*r - a*a) \
         )/ \
         (r*r - 3*M*r + 2.*a*sqrt(M*r));
}

double lnOfhTerm1(double a,
                double r, double theta, 
                double l)
{
  double Delta = computeDelta(a, r, theta);
  double Sigma = computeSigma(a, r, theta);
  double A     = computeA(a, r, theta);

  return 0.5*log( (1. + sqrt(1. + (4.*l*l*Sigma*Sigma*Delta)/ \
                                  (A*sin(theta)*A*sin(theta))
                            )
                  ) / (Sigma*Delta/A)
                );
}

double lnOfhTerm2(double a,
                double r, double theta, 
                double l)
{
  double Delta = computeDelta(a, r, theta);
  double Sigma = computeSigma(a, r, theta);
  double A     = computeA(a, r, theta);

  return -0.5*sqrt(1. + (4.*l*l*Sigma*Sigma*Delta) /
                        (A*A*sin(theta)*sin(theta))
                  );

}

double lnOfhTerm3(double a,
                double r, double theta, 
                double l)
{
  double A     = computeA(a, r, theta);
  double M = 1.;
  return -2*a*M*r*l/A;
}

double computeDelta(double a, double r, double theta)
{
  double M = 1.;
  return r*r - 2*M*r + a*a;
}

double computeSigma(double a, double r, double theta)
{
  return r*r + a*a*cos(theta)*cos(theta);
}

double computeA(double a, double r, double theta)
{
  double Delta = computeDelta(a, r, theta);

  return pow(r*r + a*a, 2.) - Delta*a*a*sin(theta)*sin(theta);
}


double computeLnOfh(double a, double r, double theta)
{
  double l = lFishboneMoncrief(a, params::PressureMaxRadius, M_PI/2.);

  double term1 = lnOfhTerm1(a, r, theta, l);
  double term2 = lnOfhTerm2(a, r, theta, l);
  double term3 = lnOfhTerm3(a, r, theta, l);

  double term1InnerEdge = lnOfhTerm1(a, params::InnerEdgeRadius, M_PI/2., l);
  double term2InnerEdge = lnOfhTerm2(a, params::InnerEdgeRadius, M_PI/2., l);
  double term3InnerEdge = lnOfhTerm3(a, params::InnerEdgeRadius, M_PI/2., l);

  return  term1 + term2 + term3 \
        - term1InnerEdge - term2InnerEdge - term3InnerEdge;

}

void timeStepper::initialConditions(int &numReads,int &numWrites)
{
  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
  //Let's ignore ArrayFire here - it's only the initial
  //conditions...
  array xCoords[3];
  geomCenter->XCoordsToxCoords(XCoords->vars, xCoords);

  const int N1g = primOld->N1Total;
  const int N2g = primOld->N2Total;
  const int N3g = primOld->N3Total;

  array& Rho = primOld->vars[vars::RHO];
  array& U = primOld->vars[vars::U];
  array& U1 = primOld->vars[vars::U1];
  array& U2 = primOld->vars[vars::U2];
  array& U3 = primOld->vars[vars::U3];
  array& B1 = primOld->vars[vars::B1];
  array& B2 = primOld->vars[vars::B2];
  array& B3 = primOld->vars[vars::B3];


  printf("Local size: %i x %i x %i\n",N1g,N2g,N3g);
  if(world_rank==0)
    printf("Running on %i procs\n",world_size);
  for(int proc=0;proc<world_size;proc++)
    {
      if(world_rank==proc)
	af_print(xCoords[directions::X1],5);
      MPI_Barrier(PETSC_COMM_WORLD);
    }
  
  double aBH = params::blackHoleSpin;
  for(int k=0;k<N3g;k++)
    for(int j=0;j<N2g;j++)
      for(int i=0;i<N1g;i++)
	{
	  const int p = i+j*N1g+k*N2g;
	  const double& r = xCoords[directions::X1](i,j,k).scalar<double>();
	  const double& theta = xCoords[directions::X2](i,j,k).scalar<double>();
	  const double& phi = xCoords[directions::X3](i,j,k).scalar<double>();
	  const double& X2 = XCoords->vars[directions::X2](i,j,k).scalar<double>();

	  const double& lapse = geomCenter->alpha(i,j,k).scalar<double>();
	  const double& beta1 = geomCenter->gCon[0][1](i,j,k).scalar<double>();
	  const double& beta2 = geomCenter->gCon[0][2](i,j,k).scalar<double>();
	  const double& beta3 = geomCenter->gCon[0][3](i,j,k).scalar<double>();


	  double lnOfh = 1.;
	  if(r>=params::InnerEdgeRadius)
	    lnOfh = computeLnOfh(aBH,r,theta);
	  
	  /* Region outside the torus */
	  if(lnOfh<0. || r<params::InnerEdgeRadius)
	    {
	      Rho(i,j,k)=params::rhoFloorInFluidElement;
	      U(i,j,k)=params::uFloorInFluidElement;
	      U1(i,j,k)=0.;
	      U2(i,j,k)=0.;
	      U3(i,j,k)=0.;
	    }
	  else
	    {
	      double h = exp(lnOfh);
	      double Gamma = params::adiabaticIndex;
	      double Kappa = params::Adiabat;

	      /* Solve for rho using the definition of h = (rho + u + P)/rho where rho
	       * here is the rest mass energy density and P = C * rho^Gamma */
	      Rho(i,j,k) = pow((h-1)*(Gamma-1.)/(Kappa*Gamma), 
			       1./(Gamma-1.));
	      U(i,j,k) =  Kappa * pow(Rho(i,j,k), Gamma)/(Gamma-1.);
	  
	      /* TODO: Should add random noise here */
	  
	      
	      /* Fishbone-Moncrief u_phi is given in the Boyer-Lindquist coordinates.
	       * Need to transform to (modified) Kerr-Schild */
	      double A = computeA(aBH, r, theta);
	      double Sigma = computeSigma(aBH, r, theta);
	      double Delta = computeDelta(aBH, r, theta);
	      double l = lFishboneMoncrief(aBH, params::PressureMaxRadius, M_PI/2.);
	      double expOfMinus2Chi = Sigma*Sigma*Delta/(A*A*sin(theta)*sin(theta)) ;
	      double uCovPhiBL = sqrt((-1. + sqrt(1. + 4*l*l*expOfMinus2Chi))/2.);
	      double uConPhiBL =   2.*aBH*r*sqrt(1. + uCovPhiBL*uCovPhiBL)
		/sqrt(A*Sigma*Delta)+ sqrt(Sigma/A)*uCovPhiBL/sin(theta);

	      double uConBL[NDIM];
	      uConBL[0] = 0.;
	      uConBL[1] = 0.;
	      uConBL[2] = 0.;
	      uConBL[3] = uConPhiBL;
	      
	      double gCovBL[NDIM][NDIM], gConBL[NDIM][NDIM];
	      double transformBLToMKS[NDIM][NDIM];
	      
	      for (int alpha=0; alpha<NDIM; alpha++)
		{
		  for (int beta=0; beta<NDIM; beta++)
		    {
		      gCovBL[alpha][beta] = 0.;
		      gConBL[alpha][beta] = 0.;
		      transformBLToMKS[alpha][beta] = 0.;
		    }
		}
	      
	      double mu = 1 + aBH*aBH*cos(theta)*cos(theta)/(r*r);
	      
	      gCovBL[0][0] = -(1. - 2./(r*mu));
	      gCovBL[0][3] = -2.*aBH*sin(theta)*sin(theta)/(r*mu);
	      gCovBL[3][0] = gCovBL[0][3];
	      gCovBL[1][1] = mu*r*r/Delta;
	      gCovBL[2][2] = r*r*mu;
	      gCovBL[3][3] = r*r*sin(theta)*sin(theta)*
		(1. + aBH*aBH/(r*r) + 2.*aBH*aBH*sin(theta)*sin(theta)/(r*r*r*mu));
	      
	      gConBL[0][0] = -1. -2.*(1 + aBH*aBH/(r*r))/(Delta*mu/r);
	      gConBL[0][3] = -2.*aBH/(r*Delta*mu);
	      gConBL[3][0] = gConBL[0][3];
	      gConBL[1][1] = Delta/(r*r*mu);
	      gConBL[2][2] = 1./(r*r*mu);
	      gConBL[3][3] = (1. - 2./(r*mu))/(sin(theta)*sin(theta)*Delta);
	      
	      transformBLToMKS[0][0] = 1.;
	      transformBLToMKS[1][1] = 1.;
	      transformBLToMKS[2][2] = 1.;
	      transformBLToMKS[3][3] = 1.;
	      transformBLToMKS[0][1] = 2.*r/Delta;
	      transformBLToMKS[3][1] = aBH/Delta; 
	      
	      /* Need to get uConBL[0] using u^mu u_mu = -1 */
	      double AA = gCovBL[0][0];
	      double BB = 2.*(gCovBL[0][1]*uConBL[1] +
			      gCovBL[0][2]*uConBL[2] +
			      gCovBL[0][3]*uConBL[3]
			      );
	      double CC = 1. + gCovBL[1][1]*uConBL[1]*uConBL[1] +
		gCovBL[2][2]*uConBL[2]*uConBL[2] +
		gCovBL[3][3]*uConBL[3]*uConBL[3] +
		2.*(gCovBL[1][2]*uConBL[1]*uConBL[2] +
		    gCovBL[1][3]*uConBL[1]*uConBL[3] +
		    gCovBL[2][3]*uConBL[2]*uConBL[3]);
	      
	      double discriminent = BB*BB - 4.*AA*CC;
	      uConBL[0] = -(BB + sqrt(discriminent))/(2.*AA);
	      
	      double uConKS[NDIM];
	      
	      for (int alpha=0; alpha<NDIM; alpha++)
		{
		  uConKS[alpha] = 0.;
		  
		  for (int beta=0; beta<NDIM; beta++)
		    {
		      uConKS[alpha] += transformBLToMKS[alpha][beta]*uConBL[beta];
		    }
		}
	      
	      /* Finally get the four-velocity in the X coordinates, which is modified
	       * Kerr-Schild */
	      double uConMKS[NDIM];
	      double rFactor = r;
	      double hFactor = M_PI + (1. - params::hSlope)*M_PI*cos(2.*M_PI*X2);
	      uConMKS[0] = uConKS[0];
	      uConMKS[1] = uConKS[1]/rFactor;
	      uConMKS[2] = uConKS[2]/hFactor;
	      uConMKS[3] = uConKS[3];
	      
	      U1(i,j,k) = uConMKS[1] + pow(lapse, 2.)*beta1*uConMKS[0];
	      U2(i,j,k) = uConMKS[2] + pow(lapse, 2.)*beta2*uConMKS[0];
	      U3(i,j,k) = uConMKS[3] + pow(lapse, 2.)*beta3*uConMKS[0];
	    }
	  B1(i,j,k) = 0.;
	  B2(i,j,k) = 0.;
	  B3(i,j,k) = 0.;
	}

  array rhoMax_af = af::max(af::max(af::max(Rho,2),1),0);
  double rhoMax = rhoMax_af.host<double>()[0];

  /* Communicate rhoMax to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  if(rhoMax < temp)
	    rhoMax = temp;
	}
      }
    else
      {
        MPI_Send(&rhoMax, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
      }
  MPI_Barrier(PETSC_COMM_WORLD);
  if (world_rank == 0)
    MPI_Bcast(&rhoMax,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  
  Rho=Rho/rhoMax;
  U=U/rhoMax;

  Rho.eval();
  U.eval();
  U1.eval();
  U2.eval();
  U3.eval();
  B1.eval();
  B2.eval();
  B3.eval();

  /* TODO : apply floors */

  /* Set magnetic field. This is MUCH easier to do using ArrayFire... */
  // Set vector potential
  const array& Rho_af = primOld->vars[vars::RHO];
  array rhoAvg = 
    (af::shift(Rho_af,1,0,0)+af::shift(Rho_af,-1,0,0)
    +af::shift(Rho_af,0,1,0)+af::shift(Rho_af,0,-1,0)
     +af::shift(Rho_af,0,0,1)+af::shift(Rho_af,0,0,-1))
    /6.;
  array zero = rhoAvg*0.;
  array Avec = af::max(rhoAvg,zero)*
    af::cos(xCoords[directions::X2]*
	    (params::MagneticLoops-1));
  Avec.eval();
 
  // Compute magnetic field 
  const array& g = geomCenter->g;
  const double dX1 = primOld->dX1;
  const double dX2 = primOld->dX2;
  primOld->vars[vars::B1] = 
    (af::shift(Avec,0,-1,0)-af::shift(Avec,0,0,0)
     +af::shift(Avec,-1,-1,0)-af::shift(Avec,-1,0,0))/
    (2.*dX2*g);
  primOld->vars[vars::B2] = 
    (af::shift(Avec,0,0,0)-af::shift(Avec,-1,0,0)
     +af::shift(Avec,0,-1,0)-af::shift(Avec,-1,-1,0))/
    (2.*dX1*g);
  primOld->vars[vars::B1].eval();
  primOld->vars[vars::B2].eval();

  // We have used ghost zones in the previous steps -> need to communicate
  // before computing global quantities
  primOld->communicate();

  // Need to set fluid element to get b^2...
  {
    elemOld->set(primOld->vars,*geomCenter,numReads,numWrites);
    const array& bSqr = elemOld->bSqr;
    const array& Pgas = elemOld->pressure;
    array PlasmaBeta = (Pgas+1.e-13)/(bSqr+1.e-18);
    array BetaMin_af = af::min(af::min(af::min(PlasmaBeta,2),1),0);
    double BFactor = BetaMin_af.host<double>()[0];
    BFactor = sqrt(BFactor/params::MinPlasmaBeta);

    /* Use MPI to find minimum over all processors */
    if (world_rank == 0) 
      {
	double temp; 
	for(int i=1;i<world_size;i++)
	  {
	    MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	    if(BFactor > temp)
	      BFactor = temp;
	  }
      }
    else
      {
        MPI_Send(&BFactor, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
      }
    MPI_Barrier(PETSC_COMM_WORLD);
    if (world_rank == 0)
      MPI_Bcast(&BFactor,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
    MPI_Barrier(PETSC_COMM_WORLD);
    
    primOld->vars[vars::B1] *= BFactor;
    primOld->vars[vars::B2] *= BFactor;
    primOld->vars[vars::B3] *= BFactor;
    primOld->vars[vars::B1].eval();
    primOld->vars[vars::B2].eval();
    primOld->vars[vars::B3].eval();
  }

  {
    elemOld->set(primOld->vars,*geomCenter,numReads,numWrites);
    const array& bSqr = elemOld->bSqr;
    const array& Pgas = elemOld->pressure;
    array PlasmaBeta = (Pgas+1.e-13)/(bSqr+1.e-18);
    array BetaMin_af = af::min(af::min(af::min(PlasmaBeta,2),1),0);
    af_print(BetaMin_af);
  }

  applyFloor(primOld,elemOld,geomCenter,numReads,numWrites);

  for (int var=0; var<vars::dof; var++) 
    {
      primIC->vars[var]=1.*primOld->vars[var];
      primIC->vars[var].eval();
    }

  af::sync();
}

void applyFloor(grid* prim, fluidElement* elem, geometry* geom, int &numReads,int &numWrites)
{
  array condition = prim->vars[vars::RHO]<1.e-5;
  prim->vars[vars::RHO] = condition*1.e-5+(1.-condition)*prim->vars[vars::RHO];

  condition = prim->vars[vars::U]<1.e-8;
  prim->vars[vars::U] = condition*1.e-8+(1.-condition)*prim->vars[vars::U];

  elem->set(prim->vars,*geom,numReads,numWrites);

  const array& bSqr = elem->bSqr;
  condition = bSqr>10.*prim->vars[vars::RHO];
  prim->vars[vars::RHO] = prim->vars[vars::RHO]*(1.-condition)+condition*0.1*bSqr;
  condition = bSqr>500.*prim->vars[vars::U];
  prim->vars[vars::U] = prim->vars[vars::U]*(1.-condition)+condition*0.002*bSqr;
  
  prim->vars[vars::RHO].eval();
  prim->vars[vars::U].eval();

  af::sync();
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{
  applyFloor(primHalfStep,elemHalfStep,geomCenter,numReads,numWrites);
  //printf("halfStepDiagnostics\n");
}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
  applyFloor(primOld,elemOld,geomCenter,numReads,numWrites);
  /* Compute the errors for the different modes */
  //printf("fullStepDiagnostics: Time = %e\n",params::Time);
  //af_print(primOld->vars[vars::RHO],8.);
}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{
};
