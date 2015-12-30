#include "physics.hpp"

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi = one;
  nu  = one;
}

fluidElement::fluidElement(const grid &prim,
                           const geometry &geom,
                           const int location
                          )
{
  /* Use this to set various fluid parameters. For ex: tau = 0.1*one etc..*/
  one = af::constant(1, 
                     prim.vars[0].dims(directions::X1),
                     prim.vars[0].dims(directions::X2),
                     prim.vars[0].dims(directions::X3),
		     f64
		     );

  set(prim, geom, location);
}

void fluidElement::set(const grid &prim,
                       const geometry &geom,
                       const int location
                      )
{
  loc = location;

  rho = prim.vars[vars::RHO] + params::rhoFloorInFluidElement;
  u   = prim.vars[vars::U  ] + params::uFloorInFluidElement;
  u1  = prim.vars[vars::U1 ];
  u2  = prim.vars[vars::U2 ];
  u3  = prim.vars[vars::U3 ];
  B1  = prim.vars[vars::B1 ];
  B2  = prim.vars[vars::B2 ];
  B3  = prim.vars[vars::B3 ];

  pressure    = (params::adiabaticIndex - 1.)*u;
  temperature = pressure/rho + params::temperatureFloorInFluidElement;

  setFluidElementParameters(geom);
  
  if (params::conduction==1)
  {
    qTilde = prim.vars[vars::Q];

    if (params::highOrderTermsConduction==1)
    {
      q = qTilde * temperature * af::sqrt(rho*chi/tau);
    }
    else
    {
      q = qTilde;
    }
  }

  if (params::viscosity==1)
  {
    deltaPTilde = prim.vars[vars::DP];

    if (params::highOrderTermsViscosity == 1)
    {
      deltaP = deltaPTilde * af::sqrt(temperature * rho * nu / tau);
    }
    else
    {
      deltaP = deltaPTilde;
    }
  }

  gammaLorentzFactor =
    af::sqrt(1 + geom.gCov[loc][1][1] * u1 * u1
               + geom.gCov[loc][2][2] * u2 * u2
               + geom.gCov[loc][3][3] * u3 * u3

             + 2*(  geom.gCov[loc][1][2] * u1 * u2
                  + geom.gCov[loc][1][3] * u1 * u3
                  + geom.gCov[loc][2][3] * u2 * u3
                 )
            );

  uCon[0] = gammaLorentzFactor/geom.alpha[loc];
  uCon[1] = u1 - gammaLorentzFactor*geom.gCon[loc][0][1]*geom.alpha[loc];
  uCon[2] = u2 - gammaLorentzFactor*geom.gCon[loc][0][2]*geom.alpha[loc];
  uCon[3] = u3 - gammaLorentzFactor*geom.gCon[loc][0][3]*geom.alpha[loc];

  for (int mu=0; mu < NDIM; mu++)
  {
    uCov[mu] =  geom.gCov[loc][mu][0] * uCon[0]
              + geom.gCov[loc][mu][1] * uCon[1]
              + geom.gCov[loc][mu][2] * uCon[2]
              + geom.gCov[loc][mu][3] * uCon[3];
  }

  bCon[0] =  B1*uCov[1] + B2*uCov[2] + B3*uCov[3];

  bCon[1] = (B1 + bCon[0] * uCon[1])/uCon[0];
  bCon[2] = (B2 + bCon[0] * uCon[2])/uCon[0];
  bCon[3] = (B3 + bCon[0] * uCon[3])/uCon[0];

  for (int mu=0; mu < NDIM; mu++)
  {
    bCov[mu] =  geom.gCov[loc][mu][0] * bCon[0]
              + geom.gCov[loc][mu][1] * bCon[1]
              + geom.gCov[loc][mu][2] * bCon[2]
              + geom.gCov[loc][mu][3] * bCon[3];
  }

  bSqr =  bCon[0]*bCov[0] + bCon[1]*bCov[1]
        + bCon[2]*bCov[2] + bCon[3]*bCov[3] + params::bSqrFloorInFluidElement;

  for (int mu : indicesToLoopOver[loc])
  {
    NUp[mu] = rho * uCon[mu];

    for (int nu=0; nu < NDIM; nu++)
    {
      TUpDown[mu][nu] =   (rho + u + pressure + bSqr)*uCon[mu]*uCov[nu]
                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)
                        - bCon[mu] * bCov[nu];

      if (params::conduction==1)
      {
        TUpDown[mu][nu] += q/af::sqrt(bSqr) * (uCon[mu]*bCov[nu] + bCon[mu]*uCov[nu]);
      }

      if (params::viscosity==1)
      {
        TUpDown[mu][nu] += (- deltaP)       
                           * (  bCon[mu] * bCov[nu]/bSqr
                              - (1./3.)*(DELTA(mu, nu) + uCon[mu]*uCov[nu])
                             );
      }
    }
  }

}

void fluidElement::computeFluxes(const geometry &geom, 
                                 const int dir,
                                 grid &flux
                                )
{
  array g = geom.g[loc];

  flux.vars[vars::RHO] = g*NUp[dir];

  flux.vars[vars::U]   = g*TUpDown[dir][0] + flux.vars[vars::RHO];

  flux.vars[vars::U1]  = g*TUpDown[dir][1];

  flux.vars[vars::U2]  = g*TUpDown[dir][2];

  flux.vars[vars::U3]  = g*TUpDown[dir][3];

  flux.vars[vars::B1]  = g*(bCon[1]*uCon[dir] - bCon[dir]*uCon[1]);

  flux.vars[vars::B2]  = g*(bCon[2]*uCon[dir] - bCon[dir]*uCon[2]);

  flux.vars[vars::B3]  = g*(bCon[3]*uCon[dir] - bCon[dir]*uCon[3]);

  if (params::conduction)
  {
    flux.vars[vars::Q] = g*(uCon[dir] * qTilde);
  }

  if (params::viscosity)
  {
    flux.vars[vars::DP] = g*(uCon[dir] * deltaPTilde);
  }

}

void fluidElement::computeSources(const geometry &geom,
				                          const fluidElement &elemOld,
                        				  const fluidElement &elemForSpatialDeriv,
                        				  const double dt,
                        				  const int useImplicitSources,
                                  grid &sources
                                 )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var] = 0.;
  }

  if (params::metric == metrics::MODIFIED_KERR_SCHILD)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int kappa=0; kappa<NDIM; kappa++)
      {
        for (int lamda=0; lamda<NDIM; lamda++)
        {
          sources.vars[vars::U + nu] +=
            geom.g[locations::CENTER]
          * TUpDown[kappa][lamda]
          * geom.gammaUpDownDown[lamda][kappa][nu];
        }
      }
    }
  }

  if (params::conduction || params::viscosity)
    {
      riemannSolver riemann(geom);

      // Non-ideal pieces
      // First, compute part of the source terms
      // shared by conduction and viscosity, i.e. 
      // u_{\mu;\nu} and u^{\mu}_{;\mu}
      array graduCov[NDIM][NDIM];     
      array divuCov = af::constant(0., rho.dims(0), rho.dims(1), rho.dims(2),f64);
      // Use this if the source terms are treated implicitly (we then
      // compute them at the beginning and end of the step)
      array graduCovOld[NDIM][NDIM];     
      array divuCovOld = af::constant(0., rho.dims(0), rho.dims(1), rho.dims(2),f64);
      for(int mu=0;mu<NDIM;mu++)
	{
	  for(int nu=0;nu<NDIM;nu++)
	    {
	      graduCov[nu][mu] = af::constant(0., rho.dims(0), rho.dims(1), rho.dims(2),f64);
	      graduCovOld[nu][mu] = af::constant(0., rho.dims(0), rho.dims(1), rho.dims(2),f64);
	    }
	  // 1) Time derivatives. Note that the old element is ghosted.
	  graduCov[0][mu] += (uCov[mu]-elemOld.uCov[mu])/dt;
	  // 2) Spatial derivatives.
	  array du = riemann.slopeMM(directions::X1,elemForSpatialDeriv.uCov[mu]);
	  graduCov[1][mu] += du;
	  du =  riemann.slopeMM(directions::X2,elemForSpatialDeriv.uCov[mu]);
          graduCov[2][mu] += du;
	  du =  riemann.slopeMM(directions::X3,elemForSpatialDeriv.uCov[mu]);
          graduCov[3][mu] += du;
	  // 3) Connection terms
	  if(useImplicitSources)
	    {
	      for(int nu=0;nu<NDIM;nu++)
		for(int mu=0;mu<NDIM;mu++)
		  {
		    graduCovOld[nu][mu]=graduCov[nu][mu];
		    for(int lambda=0;lambda<NDIM;lambda++)
		      {
			graduCov[nu][mu]+=geom.gammaUpDownDown[lambda][nu][mu]*uCov[lambda];
			graduCovOld[nu][mu]+=geom.gammaUpDownDown[lambda][nu][mu]*
			  elemOld.uCov[lambda];
		      }
		  }
	    }
	  else
	    {
	      for(int nu=0;nu<NDIM;nu++)
                for(int lambda=0;lambda<NDIM;lambda++)
                  graduCov[nu][mu]+=geom.gammaUpDownDown[lambda][nu][mu]*
		    elemForSpatialDeriv.uCov[lambda];
	    }
	}
      // 4) Compute divergence. We could compute it from the derivatives of uCon,
      //    but we already have u_{\mu;\nu}, so let's make use of it
      for(int mu=0;mu<NDIM;mu++)
	for(int nu=0;nu<NDIM;nu++)
	  {
	    divuCov += geom.gCon[mu][nu]*graduCov[mu][nu];
	    if(useImplicitSources)
	      divuCovOld += geom.gCon[mu][nu]*graduCovOld[mu][nu];
	  }
      
      // -------------------------------------
      // Now, look at viscosity-specific terms
      if(params::viscosity)
	{
	  // Compute target deltaP
	  array deltaP0 = -divuCov*rho*nu;
	  if(useImplicitSources)
	    {
	      deltaP0 = 0.5*(deltaP0 - divuCovOld*elemOld.rho*elemOld.nu);
	      for(int mu=0;mu<NDIM;mu++)
		for(int nu=0;nu<NDIM;nu++)
		  deltaP0 += 0.5*(3.*rho*nu*bCon[mu]*bCon[nu]/bSqr*graduCov[mu][nu]+
				  3.*elemOld.rho*elemOld.nu
				    *elemOld.bCon[mu]
				    *elemOld.bCon[nu]
				    /elemOld.bSqr
				    *graduCovOld[mu][nu]);
	    }
	  else
	    for(int mu=0;mu<NDIM;mu++)
	      for(int nu=0;nu<NDIM;nu++)
		deltaP0 += 3.*elemForSpatialDeriv.rho
		  *elemForSpatialDeriv.nu
		  *elemForSpatialDeriv.bCon[mu]
		  *elemForSpatialDeriv.bCon[nu]
		  /elemForSpatialDeriv.bSqr
		  *graduCov[mu][nu];
	  
	  if (params::highOrderTermsViscosity == 1)
	    deltaP0 *= af::sqrt(elemForSpatialDeriv.tau
			    /elemForSpatialDeriv.nu
			    /elemForSpatialDeriv.rho
			    /elemForSpatialDeriv.temperature);
	  
	  
	  if(useImplicitSources)
	    {
	      sources.vars[vars::DP] += (deltaP0 - 0.5*deltaPTilde 
					 - 0.5*elemOld.deltaPTilde)
		/elemForSpatialDeriv.tau;
	      if (params::highOrderTermsViscosity == 1)
		sources.vars[vars::DP] += 0.25*divuCov*deltaPTilde
		  +0.25*divuCovOld*elemOld.deltaPTilde;
	    }
	  else
	    {
	      sources.vars[vars::DP] += (deltaP0 - elemForSpatialDeriv.deltaPTilde)
		/elemForSpatialDeriv.tau;
	      if (params::highOrderTermsViscosity == 1)
		sources.vars[vars::DP] +=0.5*divuCov*elemForSpatialDeriv.deltaPTilde;
	    }
	}

      // -------------------------------------
      // Finally, look at conduction-specific terms
      if(params::conduction)
	{
	  array gradT[NDIM];
	  for(int mu=0;mu<NDIM;mu++)
	    gradT[mu] = af::constant(0., rho.dims(0), rho.dims(1), rho.dims(2),f64);
	  gradT[0] = (temperature - elemOld.temperature)/dt;
	  array dT = riemann.slopeMM(directions::X1,elemForSpatialDeriv.temperature);
	  gradT[1] = dT;
	  dT = riemann.slopeMM(directions::X2,elemForSpatialDeriv.temperature);
          gradT[2] = dT;
	  dT = riemann.slopeMM(directions::X3,elemForSpatialDeriv.temperature);
          gradT[3] = dT;
	  array q0 = af::constant(0., rho.dims(0), rho.dims(1), rho.dims(2),f64);
	  if(useImplicitSources)
	    {
	      array bnorm = af::sqrt(bSqr);
	      array bnormOld = af::sqrt(elemOld.bSqr);
	      for(int mu=0;mu<NDIM;mu++)
		{
		  q0 -= 0.5*(rho*chi*bCon[mu]/bnorm
			     +elemOld.rho
			     *elemOld.chi
			     *elemOld.bCon[mu]
			     /bnormOld)
		    *gradT[mu];
		  for(int nu=0;nu<NDIM;nu++)
		    {
		      q0 -= 0.5*(rho*chi*temperature*bCon[mu]/bnorm*uCon[nu]*graduCov[mu][nu]
				 +elemOld.rho
				 *elemOld.chi
				 *elemOld.temperature
				 *elemOld.bCon[mu]
				 /bnormOld
				 *elemOld.uCon[nu]
				 *graduCovOld[mu][nu]);
		    }
		}
	    }
	  else
	    {
	      array bnorm = af::sqrt(elemForSpatialDeriv.bSqr);
	      for(int mu=0;mu<NDIM;mu++)
		{
		  q0 -= elemForSpatialDeriv.rho
		    *elemForSpatialDeriv.chi
		    *elemForSpatialDeriv.bCon[mu]
		    /bnorm
		    *gradT[mu];
		  for(int nu=0;nu<NDIM;nu++)
		    q0 -= elemForSpatialDeriv.rho
		      *elemForSpatialDeriv.chi
		      *elemForSpatialDeriv.temperature
		      *elemForSpatialDeriv.bCon[mu]
		      /bnorm
		      *elemForSpatialDeriv.uCon[nu]
		      *graduCov[mu][nu];
		}
	    }
	  if (params::highOrderTermsConduction == 1)
	    q0 *= af::sqrt(elemForSpatialDeriv.tau
		       /elemForSpatialDeriv.chi
		       /elemForSpatialDeriv.rho)
	      /elemForSpatialDeriv.temperature;
	  
	  
	  if(useImplicitSources)
	    {
	      sources.vars[vars::Q] += (q0 - 0.5*qTilde 
					 - 0.5*elemOld.qTilde)
		/elemForSpatialDeriv.tau;
	      if (params::highOrderTermsConduction == 1)
		sources.vars[vars::Q] += 0.25*divuCov*qTilde
		  +0.25*divuCovOld*elemOld.qTilde;
	    }
	  else
	    {
	      sources.vars[vars::Q] += (q0 - elemForSpatialDeriv.qTilde)
		/elemForSpatialDeriv.tau;
	      if (params::highOrderTermsConduction == 1)
		sources.vars[vars::Q] +=0.5*divuCov*elemForSpatialDeriv.qTilde;
	    }

	}
    }
}
