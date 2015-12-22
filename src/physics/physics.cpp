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
  set(prim, geom, location);
}

void fluidElement::set(const grid &prim,
                       const geometry &geom,
                       const int location
                      )
{
  /* Set locations along axisymmetric direction to that of the center of the
   * grid zone */
  if (location==locations::FRONT || 
      location==locations::BACK
     )
  {
    loc = locations::CENTER;
  }
  else
  {
    loc = location;
  }

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

  one = af::constant(1, rho.dims(0), rho.dims(1), rho.dims(2));
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

  uCon[0] = gamma/geom.alpha[loc];
  uCon[1] = u1 - gamma*geom.gCon[loc][0][1]*geom.alpha[loc];
  uCon[2] = u2 - gamma*geom.gCon[loc][0][2]*geom.alpha[loc];
  uCon[3] = u3 - gamma*geom.gCon[loc][0][3]*geom.alpha[loc];

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
        TUpDown[mu][nu] += q/sqrt(bSqr) * (uCon[mu]*bCov[nu] + bCon[mu]*uCov[nu]);
      }

      if (params::viscosity==1)
      {
        TUpDown[mu][nu] += - deltaP       
                           * (  bCon[mu] * bCov[nu]/bSqr
                              + (1./3.)*(DELTA(mu, nu) + uCon[mu]*uCov[nu])
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

}
