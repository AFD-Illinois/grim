#include "physics.hpp"

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
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
	dtuCov[mu] = zero;
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

  soundSpeed  = af::sqrt( params::adiabaticIndex*pressure
                         /(rho+params::adiabaticIndex*u)
                        );

  setFluidElementParameters(geom);
  
  if (params::conduction==1)
  {
    qTilde = prim.vars[vars::Q];

    if (params::highOrderTermsConduction==1)
    {
      q = qTilde * temperature * af::sqrt(rho*chi_emhd/tau);
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
      deltaP = deltaPTilde * af::sqrt(temperature * rho * nu_emhd / tau);
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

void fluidElement::computeMinMaxCharSpeeds(const geometry &geom,
			     const int dir,
			     array &MinSpeed,
			     array &MaxSpeed)
{
  array zero = 0.*one;
  
  // Alven speed
  array cAlvenSqr = bSqr/(rho+params::adiabaticIndex*u+bSqr);
  // Hydro sound speed
  array csSqr = soundSpeed*soundSpeed;
  // Approximate contribution from dP
  array cVisSqr = zero;
  if(params::viscosity)
    cVisSqr = 4./3./(rho+params::adiabaticIndex*u)*rho*nu_emhd/tau;
  // Approximate contribution from Q
  array cConSqr = zero;
  if(params::conduction)
    cConSqr = (params::adiabaticIndex-1.)*chi_emhd/tau;
  cConSqr = 0.5*(csSqr+cConSqr+af::sqrt(csSqr*csSqr+cConSqr*cConSqr));
  csSqr = cConSqr + cVisSqr;

  // Combine speeds using relativistic velocity additions
  csSqr = csSqr + cAlvenSqr - csSqr*cAlvenSqr; 
  
  array condition = csSqr > one;
  csSqr = csSqr*(one-condition)+condition;
  
  int sdir = 0;
  switch(dir)
    {
    case directions::X1:
      sdir=1; break;
    case directions::X2:
      sdir=2; break;
    case directions::X3:
      sdir=3; break;
    }
  
  array ACov[NDIM], ACon[NDIM];
  array BCov[NDIM], BCon[NDIM];
  for (int mu=0; mu<NDIM; mu++)
    {
      ACov[mu] = zero;
      BCov[mu] = zero;
    }
  ACov[sdir] = one;
  BCov[0] = one;
  for (int mu=0; mu<NDIM; mu++)
    {
      ACon[mu]=zero;
      BCon[mu]=zero;
      for(int nu=0;nu<NDIM; nu++)
	{
	  ACon[mu]+=geom.gCon[mu][nu]*ACov[nu];
	  BCon[mu]+=geom.gCon[mu][nu]*BCov[nu];
	}
    }
  array ASqr = zero;
  array BSqr = zero;
  array ADotU = zero;
  array BDotU = zero;
  array ADotB = zero;
  for (int mu=0; mu<NDIM; mu++)
    {
      ASqr += ACov[mu]*ACon[mu];
      BSqr += BCov[mu]*BCon[mu];
      ADotU += ACov[mu]*uCon[mu];
      BDotU += BCov[mu]*uCon[mu];
      ADotB += ACov[mu]*BCon[mu];
    }
  
  array A = (BDotU*BDotU) - (BSqr + BDotU*BDotU)*csSqr;
  array B = 2.*(ADotU*BDotU - (ADotB + ADotU*BDotU)*csSqr);
  array C = ADotU*ADotU - (ASqr + ADotU*ADotU)*csSqr;
  array discr = af::sqrt(B*B - 4.*A*C);
  
  MinSpeed = -(-B + discr)/2./A;
  MaxSpeed = -(-B - discr)/2./A;
  
  condition = (MinSpeed > -1.e-15);
  MinSpeed = MinSpeed*(1-condition)-1.e-15*condition;

  condition = (MaxSpeed < 1.e-15);
  MaxSpeed = MaxSpeed*(1-condition)+1.e-15*condition;
}

void fluidElement::computeTimeDerivSources(const geometry &geom,
                                          const fluidElement &elemOld,
                                          const fluidElement &elemNew,
                                          const double dt,
                                          grid &sources)
{
  for (int var=0; var<vars::dof; var++)
    {
      sources.vars[var] = 0.;
    }
  if (params::conduction || params::viscosity)
  {
    // Non-ideal pieces. Note that we only compute
    // the terms proportional to time derivatives. 

    // First, compute part of the source terms
    // shared by conduction and viscosity, i.e. 
    // u_{\mu;\nu} and u^{\mu}_{;\mu}
    // We only include the time derivative here!
    for(int mu=0;mu<NDIM;mu++)
      dtuCov[mu] = (elemNew.uCov[mu]-elemOld.uCov[mu])/dt;
    // Compute divergence. We could compute it from the derivatives of uCon,
    //    but we already have u_{\mu;\nu}, so let's make use of it
    // Naturally, this is not truly the divergence. Only the part proportional
    // to dt_u
    divuCov = 0.;
    for(int nu=0;nu<NDIM;nu++)
      divuCov += geom.gCon[0][nu]*(dtuCov[nu]);
      
    // -------------------------------------
    // Now, look at viscosity-specific terms
    if(params::viscosity)
      {
	// Compute target deltaP (time deriv part)
	deltaP0 = -divuCov*rho*nu_emhd;	    
	for(int nu=0;nu<NDIM;nu++)
	  deltaP0 += 3. * rho*nu_emhd*bCon[0]*bCon[nu]
	    /bSqr*dtuCov[nu];
	if (params::highOrderTermsViscosity == 1)
	  deltaP0 *= af::sqrt(tau/rho/nu_emhd/temperature);
	
	//Note on sign: we put the sources on the LHS when
	//computing the residual!
	sources.vars[vars::DP] -= geom.g[locations::CENTER]*(deltaP0)/tau;

	if (params::highOrderTermsViscosity == 1)
	  sources.vars[vars::DP] -= 0.5*geom.g[locations::CENTER]*divuCov*deltaPTilde;
      } /* End of viscosity specific terms */
    
    // -------------------------------------
    // Finally, look at conduction-specific terms (time deriv terms)
    if(params::conduction)
      {
	q0=0.;
	q0 -=  rho*chi_emhd*bCon[0]
	  /bNorm*(elemNew.temperature - elemOld.temperature)/dt;
	for(int nu=0;nu<NDIM;nu++)
	  q0 -=  rho*chi_emhd*temperature*
	    bCon[nu]/bNorm*uCon[0]*dtuCov[nu];
	
	if (params::highOrderTermsConduction == 1)
	  q0 *=  af::sqrt(tau/rho/chi_emhd)/temperature; 
	    
	//Note on sign: we put the sources on the LHS when
	//computing the residual!
	sources.vars[vars::Q] -= geom.g[locations::CENTER]*(q0)/tau;
	
	if (params::highOrderTermsConduction == 1)
	  sources.vars[vars::Q] -= 0.5*geom.g[locations::CENTER]*divuCov*qTilde;
      } /* End of conduction */
  } /* End of EMHD: viscosity || conduction */
}


void fluidElement::computeImplicitSources(const geometry &geom,
					  grid &sources)
{
  for (int var=0; var<vars::dof; var++)
    {
      sources.vars[var] = 0.;
    }
  if (params::conduction || params::viscosity)
  {
    // Non-ideal pieces. Note that we only compute
    // the terms treated implicitly. 
    // Look at viscosity-specific terms
    if(params::viscosity)
      {
	//Note on sign: we put the sources on the LHS when
	//computing the residual!
	sources.vars[vars::DP] -= geom.g[locations::CENTER]*(- deltaPTilde)/tau;
      } /* End of viscosity specific terms */
    
    // -------------------------------------
    // Finally, look at conduction-specific terms (implicit terms)
    if(params::conduction)
      {
	//Note on sign: we put the sources on the LHS when
	//computing the residual!
	sources.vars[vars::Q] -= geom.g[locations::CENTER]*(- qTilde)/tau;
      } /* End of conduction */
  } /* End of EMHD: viscosity || conduction */
  
}

void fluidElement::computeExplicitSources(const geometry &geom,
					  grid &sources
					  )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var] = 0.;
  }

  //Note on sign: residual computation places
  // the source terms on the LHS of the equation!
  // All ideal MHD terms are treated explicitly.
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
    // Non-ideal pieces (explicit parts)
    // First, compute part of the source terms
    // shared by conduction and viscosity, i.e. 
    // u_{\mu;\nu} and u^{\mu}_{;\mu}
    // Note that the derivatives are all precomputed 
    // in computeEMHDGradients

    // Compute divergence. We could compute it from the derivatives of uCon,
    //    but we already have u_{\mu;\nu}, so let's make use of it
    // Note that this does NOT include terms proportional to dt_u
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
	// Compute target deltaP (explicit part)
	deltaP0 = -divuCov*rho*nu_emhd;
	
	for(int mu=0;mu<NDIM;mu++)
	  for(int nu=0;nu<NDIM;nu++)
	    deltaP0 += 3. * rho * nu_emhd * bCon[mu] * bCon[nu]
	      /bSqr*graduCov[mu][nu];
	if (params::highOrderTermsViscosity == 1)
	  deltaP0 *= af::sqrt(tau/rho/nu_emhd/temperature);
	
	//Note on sign: we put the sources on the LHS when
	//computing the residual!
	// The damping term proportional to deltaPTilde is in the implicit sector.
	sources.vars[vars::DP] -=   geom.g[locations::CENTER]*(deltaP0)/tau;
	if (params::highOrderTermsViscosity == 1)
	  {
	    sources.vars[vars::DP] -= 0.5*geom.g[locations::CENTER]*divuCov*deltaPTilde;
	  }
      } /* End of viscosity specific terms */
    
    // -------------------------------------
    // Finally, look at conduction-specific terms (explicit part)
    if(params::conduction)
      {
	q0 = 0.;
	//q0 is not exactly targetQ, as the time derivative parts
	// are in the implicit sector
	for(int mu=0;mu<NDIM;mu++)
	  {
	    q0 -= rho*chi_emhd*bCon[mu]/bNorm*gradT[mu];		
	    for(int nu=0;nu<NDIM;nu++)
	      q0 -=   rho*chi_emhd*temperature*bCon[nu]/bNorm*uCon[mu]*graduCov[mu][nu];
	  }
	if (params::highOrderTermsConduction == 1)
	  q0 *=  af::sqrt(  tau/rho/chi_emhd)/temperature;
	//Note on sign: we put the sources on the LHS when
	//computing the residual!
	// The damping term proportional to qTilde is in the implicit sector. 
	sources.vars[vars::Q] -= geom.g[locations::CENTER]*(q0)/tau;
	if (params::highOrderTermsConduction == 1)
	  sources.vars[vars::Q] -= 0.5*geom.g[locations::CENTER]*divuCov*qTilde;
      } /* End of conduction */
  } /* End of EMHD: viscosity || conduction */
}

void fluidElement::computeEMHDGradients(const geometry &geom)
{
  if (params::conduction || params::viscosity)
  {
    double dX1 = geom.XCoordsGrid->dX1;
    double dX2 = geom.XCoordsGrid->dX2;
    double dX3 = geom.XCoordsGrid->dX3;
    
    for(int mu=0;mu<NDIM;mu++)
	  {
	    //Time derivative needs to be reset for reach residual computation,
	    // so not computed here.
	    graduCov[0][mu] = 0.;
	  
	    array du = reconstruction::slope(directions::X1,dX1,uCov[mu]);
	    graduCov[1][mu] = du;

	    if(params::dim>1)
	    {
	      du =  reconstruction::slope(directions::X2,dX2,uCov[mu]);
	      graduCov[2][mu] = du;
	    }
	  
	    if(params::dim>2)
	    {
	      du =  reconstruction::slope(directions::X3,dX3,uCov[mu]);
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

	    array dT = reconstruction::slope(directions::X1,dX1,temperature);
	    gradT[1] = dT;

	    if(params::dim>1)
	    {
	      dT = reconstruction::slope(directions::X2,dX2,temperature);
	      gradT[2] = dT;
	    }

	    if(params::dim>2)
	    {
	      dT = reconstruction::slope(directions::X3,dX3,temperature);
	      gradT[3] = dT;
	    }
	  } /* End of conduction specific terms */

  }
}
