#include "physics.hpp"

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi = soundSpeed*soundSpeed*tau;
  nu  = soundSpeed*soundSpeed*tau;
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

  array zero = 0.*one;

  /* Allocate memory for gradients used in EMHD */
  if (params::conduction || params::viscosity)
  {
    divuCov = zero;
    
    for(int mu=0;mu<NDIM;mu++)
	  {
	    gradT[mu] = zero;
	    
      for(int nu=0;nu<NDIM;nu++)
      {
	      graduCov[nu][mu] = zero;
      }
	  }
      
    deltaP0 = zero;
    q0 = zero;
  }

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

  soundSpeed = af::sqrt(params::adiabaticIndex*pressure/(rho+params::adiabaticIndex*u));
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
  bNorm = af::sqrt(bSqr);

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
        TUpDown[mu][nu] += q/bNorm * (uCon[mu]*bCov[nu] + bCon[mu]*uCov[nu]);
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

  //Note on sign: residual computation places
  // the source terms on the LHS of the equation!
  if (params::metric != metrics::MINKOWSKI)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int kappa=0; kappa<NDIM; kappa++)
      {
        for (int lamda=0; lamda<NDIM; lamda++)
        {
          sources.vars[vars::U + nu] -=
            geom.g[locations::CENTER]
          * TUpDown[kappa][lamda]
          * geom.gammaUpDownDown[lamda][kappa][nu];
        }
      }
    }
  }

  if (params::conduction || params::viscosity)
  {
    // Non-ideal pieces
    // First, compute part of the source terms
    // shared by conduction and viscosity, i.e. 
    // u_{\mu;\nu} and u^{\mu}_{;\mu}
    for(int mu=0;mu<NDIM;mu++)
    {  
      for(int nu=0;nu<NDIM;nu++)
      {
        graduCov[mu][nu] = elemForSpatialDeriv.graduCov[mu][nu];
      }
    }
    // Add time derivatives, which were not precomputed of course... 
    for(int mu=0;mu<NDIM;mu++)
    {
      graduCov[0][mu] += (uCov[mu]-elemOld.uCov[mu])/dt;
    }
      
    // Compute divergence. We could compute it from the derivatives of uCon,
    //    but we already have u_{\mu;\nu}, so let's make use of it
    divuCov = 0.;
    for(int mu=0;mu<NDIM;mu++)
    {  
      for(int nu=0;nu<NDIM;nu++)
      {
	      divuCov += geom.gCon[mu][nu]*graduCov[mu][nu];
      }
    }
      
    // -------------------------------------
    // Now, look at viscosity-specific terms
    if(params::viscosity)
	  {
	    // Compute target deltaP
	    deltaP0 = -divuCov*elemForSpatialDeriv.rho
	      *elemForSpatialDeriv.nu;
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		for(int nu=0;nu<NDIM;nu++)
		  {
		    deltaP0 += 3. * elemForSpatialDeriv.rho
		      * elemForSpatialDeriv.nu
		      * elemForSpatialDeriv.bCon[mu]
		      * elemForSpatialDeriv.bCon[nu]
		      / elemForSpatialDeriv.bSqr
		      * graduCov[mu][nu];
		  }
	      }
	    
	    if (params::highOrderTermsViscosity == 1)
	      {
		deltaP0 *= af::sqrt(  elemForSpatialDeriv.tau
				      / elemForSpatialDeriv.nu
				      / elemForSpatialDeriv.rho
				      /elemForSpatialDeriv.temperature
				      );
	      }
	    
	    //Note on sign: we put the sources on the LHS when
	    //computing the residual!
	    if(useImplicitSources)
	      {
		sources.vars[vars::DP] -= geom.g[locations::CENTER]*(deltaP0 - 0.5*(deltaPTilde + elemOld.deltaPTilde)
					   )
		  / elemForSpatialDeriv.tau;
		
		if (params::highOrderTermsViscosity == 1)
		  {
		    sources.vars[vars::DP] -= 0.25*geom.g[locations::CENTER]*divuCov*(deltaPTilde+elemOld.deltaPTilde);
		  }
	      }
	    else
	      {
		sources.vars[vars::DP] -=   geom.g[locations::CENTER]*(deltaP0 - elemForSpatialDeriv.deltaPTilde)
		  / elemForSpatialDeriv.tau;
		
		if (params::highOrderTermsViscosity == 1)
		  {
		    sources.vars[vars::DP] -= 0.5*geom.g[locations::CENTER]*divuCov*elemForSpatialDeriv.deltaPTilde;
		  }
	      }
	    
	  } /* End of viscosity specific terms */
    
    // -------------------------------------
    // Finally, look at conduction-specific terms
    if(params::conduction)
      {
	    gradT[0] = (temperature - elemOld.temperature)/dt;
	    
	    for(int mu=1;mu<NDIM;mu++)
	      {
		gradT[mu] = elemForSpatialDeriv.gradT[mu];
	      }
	    
	    q0 = 0.;
	    for(int mu=0;mu<NDIM;mu++)
	      {
		q0 -=   elemForSpatialDeriv.rho
		  * elemForSpatialDeriv.chi
		  * elemForSpatialDeriv.bCon[mu]
		  / elemForSpatialDeriv.bNorm
		  * gradT[mu];
		
		for(int nu=0;nu<NDIM;nu++)
		  {
		    q0 -=   elemForSpatialDeriv.rho
          		  * elemForSpatialDeriv.chi
		      * elemForSpatialDeriv.temperature
		      * elemForSpatialDeriv.bCon[mu]
		      / elemForSpatialDeriv.bNorm
		      * elemForSpatialDeriv.uCon[nu]
		      * graduCov[mu][nu];
		  }
	      }
	    if (params::highOrderTermsConduction == 1)
	      {
		q0 *=  af::sqrt(  elemForSpatialDeriv.tau
				  / elemForSpatialDeriv.chi
				  / elemForSpatialDeriv.rho
				  )
		  / elemForSpatialDeriv.temperature;
	      }
	    
	    //Note on sign: we put the sources on the LHS when
	    //computing the residual!
	    if(useImplicitSources)
	      {
		sources.vars[vars::Q] -= geom.g[locations::CENTER]*(q0 - 0.5*(qTilde + elemOld.qTilde)
					  )
		  / elemForSpatialDeriv.tau;
		
		if (params::highOrderTermsConduction == 1)
		  {
		    sources.vars[vars::Q] -= 0.25*geom.g[locations::CENTER]*divuCov*(qTilde+elemOld.qTilde);
		  }
	      }
	    else
	      {
		sources.vars[vars::Q] -=   geom.g[locations::CENTER]*(q0 - elemForSpatialDeriv.qTilde)
		  / elemForSpatialDeriv.tau;
		if (params::highOrderTermsConduction == 1)
		  {
		    sources.vars[vars::Q] -= 0.5*geom.g[locations::CENTER]*divuCov*elemForSpatialDeriv.qTilde;
		  }
	      }
	    
      } /* End of conduction */
    
  } /* End of EMHD: viscosity || conduction */
}

void fluidElement::computeEMHDGradients(const geometry &geom)
{
  if (params::conduction || params::viscosity)
  {
    double dX1 = geom.xCoordsGrid->dX1;
    double dX2 = geom.xCoordsGrid->dX2;
    double dX3 = geom.xCoordsGrid->dX3;
    
    for(int mu=0;mu<NDIM;mu++)
	  {
	    //Time derivative needs to be reset for reach residual computation,
	    // so not computed here.
	    graduCov[0][mu] = 0.;
	  
	    array du = reconstruction::slopeMM(directions::X1,dX1,uCov[mu]);
	    graduCov[1][mu] = du;

	    if(params::dim>1)
	    {
	      du =  reconstruction::slopeMM(directions::X2,dX2,uCov[mu]);
	      graduCov[2][mu] = du;
	    }
	  
	    if(params::dim>2)
	    {
	      du =  reconstruction::slopeMM(directions::X3,dX3,uCov[mu]);
	      graduCov[3][mu] = du;
	    }

	    for(int nu=0;nu<NDIM;nu++)
	      {  
		for(int lambda=0;lambda<NDIM;lambda++)
		  {
		    graduCov[nu][mu]+=geom.gammaUpDownDown[lambda][nu][mu]*uCov[lambda];
		  }
	      }
	  }
    
    if(params::conduction)
	  {
	    /* Time derivative not computed here */
	    gradT[0] = 0.;

	    array dT = reconstruction::slopeMM(directions::X1,dX1,temperature);
	    gradT[1] = dT;

	    if(params::dim>1)
	    {
	      dT = reconstruction::slopeMM(directions::X2,dX2,temperature);
	      gradT[2] = dT;
	    }

	    if(params::dim>2)
	    {
	      dT = reconstruction::slopeMM(directions::X3,dX3,temperature);
	      gradT[3] = dT;
	    }
	  } /* End of conduction specific terms */

  }
}
