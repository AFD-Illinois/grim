#include "geometry.hpp"

geometry::geometry()
{
  numGhost  = params::numGhost;

  xCoordsGrid     = new grid(LOCATIONS*NDIM);
  XCoordsGrid     = new grid(LOCATIONS*NDIM);

  alphaGrid       = new grid(AXISYM_LOCATIONS);
  gGrid           = new grid(AXISYM_LOCATIONS);
  gCovGrid        = new grid(AXISYM_LOCATIONS*NDIM*NDIM);
  gConGrid        = new grid(AXISYM_LOCATIONS*NDIM*NDIM);
  connectionGrid  = new grid(NDIM*NDIM*NDIM);

  /* Indices go from [-numGhost, N + numGhost) in each direction */
  array indices = 
    af::range(af::dim4(xCoordsGrid->N1Local + 2*numGhost,
                       xCoordsGrid->N2Local + 2*numGhost,
                       xCoordsGrid->N3Local + 2*numGhost)
             ) - numGhost;

  XCoordsGrid->vars[0 + NDIM*locations::CENTER] = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::RIGHT]  = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::LEFT]   = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::BOTTOM] = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::TOP]    = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::FRONT]  = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::BACK]   = 0.;

  /* X1 coordinate at all the locations */
  XCoordsGrid->vars[1 + NDIM*locations::CENTER] = 
    params::X1Start + (indices + 0.5)*gridParams::dX1;

  XCoordsGrid->vars[1 + NDIM*locations::LEFT]   = 
    params::X1Start + (indices      )*gridParams::dX1;

  XCoordsGrid->vars[1 + NDIM*locations::RIGHT]  = 
    params::X1Start + (indices + 1. )*gridParams::dX1;

  XCoordsGrid->vars[1 + NDIM*locations::BOTTOM] = 
    XCoordsGrid->vars[1 + NDIM*locations::CENTER];

  XCoordsGrid->vars[1 + NDIM*locations::TOP]    = 
    XCoordsGrid->vars[1 + NDIM*locations::CENTER];

  XCoordsGrid->vars[1 + NDIM*locations::FRONT]  = 
    XCoordsGrid->vars[1 + NDIM*locations::CENTER];

  XCoordsGrid->vars[1 + NDIM*locations::BACK]   = 
    XCoordsGrid->vars[1 + NDIM*locations::CENTER];

  /* X2 coordinate at all the locations */
  XCoordsGrid->vars[2 + NDIM*locations::CENTER] = 
    params::X2Start + (indices + 0.5)*gridParams::dX2;

  XCoordsGrid->vars[2 + NDIM*locations::BOTTOM] = 
    params::X2Start + (indices      )*gridParams::dX2;

  XCoordsGrid->vars[2 + NDIM*locations::TOP]    = 
    params::X2Start + (indices + 1. )*gridParams::dX2;

  XCoordsGrid->vars[2 + NDIM*locations::LEFT]   = 
    XCoordsGrid->vars[2 + NDIM*locations::CENTER];

  XCoordsGrid->vars[2 + NDIM*locations::RIGHT]  = 
    XCoordsGrid->vars[2 + NDIM*locations::CENTER];

  XCoordsGrid->vars[2 + NDIM*locations::FRONT]  = 
    XCoordsGrid->vars[2 + NDIM*locations::CENTER];

  XCoordsGrid->vars[2 + NDIM*locations::BACK]   = 
    XCoordsGrid->vars[2 + NDIM*locations::CENTER];

  /* X3 coordinate at all the locations */
  XCoordsGrid->vars[3 + NDIM*locations::CENTER] = 
    params::X3Start + (indices + 0.5)*gridParams::dX3;

  XCoordsGrid->vars[3 + NDIM*locations::BACK]   = 
    params::X3Start + (indices      )*gridParams::dX3;

  XCoordsGrid->vars[3 + NDIM*locations::FRONT]  = 
    params::X3Start + (indices + 1. )*gridParams::dX3;

  XCoordsGrid->vars[3 + NDIM*locations::LEFT]   = 
    XCoordsGrid->vars[3 + NDIM*locations::CENTER];

  XCoordsGrid->vars[3 + NDIM*locations::RIGHT]  = 
    XCoordsGrid->vars[3 + NDIM*locations::CENTER];

  XCoordsGrid->vars[3 + NDIM*locations::BOTTOM] = 
    XCoordsGrid->vars[3 + NDIM*locations::CENTER];

  XCoordsGrid->vars[3 + NDIM*locations::TOP]    = 
    XCoordsGrid->vars[3 + NDIM*locations::CENTER];

  for (int loc: allLocations)
  {
    for (int mu=0; mu<NDIM; mu++)
    {
      XCoords[loc][mu] = XCoordsGrid->vars[mu + NDIM*loc];
    }

    XCoordsToxCoords(XCoords[loc], xCoords[loc]);

    for (int mu=0; mu<NDIM; mu++)
    {
      xCoordsGrid->vars[mu + NDIM*loc] = xCoords[loc][mu];
    }
  }

  /* Use the following array to allocate other arrays */
  zero = 0.*XCoords[locations::CENTER][0];

  for (int loc : axiSymmetricLocations)
  {
    /* Allocate space */
    g[loc]     = zero;
    array gDet = zero;
    for (int mu=0; mu<NDIM; mu++)
    {
      for (int nu=0; nu<NDIM; nu++)
      {
        gCov[loc][mu][nu] = zero;
        gCon[loc][mu][nu] = zero;
      }
    }
    
    /* Set gCov[loc] */
    setgCovInXCoords(XCoords[loc], gCov[loc]);

    setgDetAndgConFromgCov(gCov[loc], gDet, gCon[loc]);
    g[loc] = af::sqrt(-gDet);

    alpha[loc] = 1./af::sqrt(-gCon[loc][0][0]);

    gGrid->vars[loc]     = g[loc];
    alphaGrid->vars[loc] = alpha[loc];
    for (int mu=0; mu<NDIM; mu++)
    {
      for (int nu=0; nu<NDIM; nu++)
      {
        gCovGrid->vars[mu + NDIM*loc] = gCov[loc][mu][nu];
        gConGrid->vars[mu + NDIM*loc] = gCon[loc][mu][nu];
      }
    }
  }

}

void geometry::setgDetAndgConFromgCov(const array gCov[NDIM][NDIM],
                                      array &gDet,
                                      array gCon[NDIM][NDIM]
                                     )
{
  gDet = 
        gCov[0][0]*gCov[1][1]*gCov[2][2]*gCov[3][3] 
      + gCov[0][0]*gCov[1][2]*gCov[2][3]*gCov[3][1]
      + gCov[0][0]*gCov[1][3]*gCov[2][1]*gCov[3][2] 
      + gCov[0][1]*gCov[1][0]*gCov[2][3]*gCov[3][2] 
      + gCov[0][1]*gCov[1][2]*gCov[2][0]*gCov[3][3] 
      + gCov[0][1]*gCov[1][3]*gCov[2][2]*gCov[3][0] 
      + gCov[0][2]*gCov[1][0]*gCov[2][1]*gCov[3][3] 
      + gCov[0][2]*gCov[1][1]*gCov[2][3]*gCov[3][0] 
      + gCov[0][2]*gCov[1][3]*gCov[2][0]*gCov[3][1] 
      + gCov[0][3]*gCov[1][0]*gCov[2][2]*gCov[3][1]
      + gCov[0][3]*gCov[1][1]*gCov[2][0]*gCov[3][2] 
      + gCov[0][3]*gCov[1][2]*gCov[2][1]*gCov[3][0]
      - gCov[0][0]*gCov[1][1]*gCov[2][3]*gCov[3][2]
      - gCov[0][0]*gCov[1][2]*gCov[2][1]*gCov[3][3]
      - gCov[0][0]*gCov[1][3]*gCov[2][2]*gCov[3][1]
      - gCov[0][1]*gCov[1][0]*gCov[2][2]*gCov[3][3]
      - gCov[0][1]*gCov[1][2]*gCov[2][3]*gCov[3][0]
      - gCov[0][1]*gCov[1][3]*gCov[2][0]*gCov[3][2]
      - gCov[0][2]*gCov[1][0]*gCov[2][3]*gCov[3][1]
      - gCov[0][2]*gCov[1][1]*gCov[2][0]*gCov[3][3]
      - gCov[0][2]*gCov[1][3]*gCov[2][1]*gCov[3][0]
      - gCov[0][3]*gCov[1][0]*gCov[2][1]*gCov[3][2]
      - gCov[0][3]*gCov[1][1]*gCov[2][2]*gCov[3][0]
      - gCov[0][3]*gCov[1][2]*gCov[2][0]*gCov[3][1];

  gCon[0][0] = 
      (  gCov[1][1]*gCov[2][2]*gCov[3][3]
       + gCov[1][2]*gCov[2][3]*gCov[3][1]
       + gCov[1][3]*gCov[2][1]*gCov[3][2]
       - gCov[1][1]*gCov[2][3]*gCov[3][2]
       - gCov[1][2]*gCov[2][1]*gCov[3][3]
       - gCov[1][3]*gCov[2][2]*gCov[3][1])/gDet;

  gCon[0][1] = 
      (  gCov[0][1]*gCov[2][3]*gCov[3][2]
       + gCov[0][2]*gCov[2][1]*gCov[3][3]
       + gCov[0][3]*gCov[2][2]*gCov[3][1]
       - gCov[0][1]*gCov[2][2]*gCov[3][3]
       - gCov[0][2]*gCov[2][3]*gCov[3][1] 
       - gCov[0][3]*gCov[2][1]*gCov[3][2])/gDet;

  gCon[0][2] = 
      (  gCov[0][1]*gCov[1][2]*gCov[3][3]
       + gCov[0][2]*gCov[1][3]*gCov[3][1]
       + gCov[0][3]*gCov[1][1]*gCov[3][2]
       - gCov[0][1]*gCov[1][3]*gCov[3][2]
       - gCov[0][2]*gCov[1][1]*gCov[3][3]
       - gCov[0][3]*gCov[1][2]*gCov[3][1])/gDet;

  gCon[0][3] = 
      (  gCov[0][1]*gCov[1][3]*gCov[2][2]
       + gCov[0][2]*gCov[1][1]*gCov[2][3]
       + gCov[0][3]*gCov[1][2]*gCov[2][1]
       - gCov[0][1]*gCov[1][2]*gCov[2][3]
       - gCov[0][2]*gCov[1][3]*gCov[2][1]
       - gCov[0][3]*gCov[1][1]*gCov[2][2])/gDet;

  gCon[1][0] = gCon[0][1];
  
  gCon[1][1] = 
      (  gCov[0][0]*gCov[2][2]*gCov[3][3]
       + gCov[0][2]*gCov[2][3]*gCov[3][0]
       + gCov[0][3]*gCov[2][0]*gCov[3][2]
       - gCov[0][0]*gCov[2][3]*gCov[3][2]
       - gCov[0][2]*gCov[2][0]*gCov[3][3]
       - gCov[0][3]*gCov[2][2]*gCov[3][0])/gDet;

  gCon[1][2] = 
      (  gCov[0][0]*gCov[1][3]*gCov[3][2]
       + gCov[0][2]*gCov[1][0]*gCov[3][3]
       + gCov[0][3]*gCov[1][2]*gCov[3][0]
       - gCov[0][0]*gCov[1][2]*gCov[3][3]
       - gCov[0][2]*gCov[1][3]*gCov[3][0]
       - gCov[0][3]*gCov[1][0]*gCov[3][2])/gDet;

  gCon[1][3] = 
      (  gCov[0][0]*gCov[1][2]*gCov[2][3]
       + gCov[0][2]*gCov[1][3]*gCov[2][0]
       + gCov[0][3]*gCov[1][0]*gCov[2][2]
       - gCov[0][0]*gCov[1][3]*gCov[2][2]
       - gCov[0][2]*gCov[1][0]*gCov[2][3]
       - gCov[0][3]*gCov[1][2]*gCov[2][0])/gDet;

  gCon[2][0] = gCon[0][2];
  gCon[2][1] = gCon[1][2];

  gCon[2][2] =
      (  gCov[0][0]*gCov[1][1]*gCov[3][3]
       + gCov[0][1]*gCov[1][3]*gCov[3][0]
       + gCov[0][3]*gCov[1][0]*gCov[3][1]
       - gCov[0][0]*gCov[1][3]*gCov[3][1]
       - gCov[0][1]*gCov[1][0]*gCov[3][3]
       - gCov[0][3]*gCov[1][1]*gCov[3][0])/gDet;

  gCon[2][3] =
      (  gCov[0][0]*gCov[1][3]*gCov[2][1]
       + gCov[0][1]*gCov[1][0]*gCov[2][3]
       + gCov[0][3]*gCov[1][1]*gCov[2][0]
       - gCov[0][0]*gCov[1][1]*gCov[2][3]
       - gCov[0][1]*gCov[1][3]*gCov[2][0]
       - gCov[0][3]*gCov[1][0]*gCov[2][1])/gDet;

  gCon[3][0] = gCon[0][3];
  gCon[3][1] = gCon[1][3];
  gCon[3][2] = gCon[2][3];

  gCon[3][3] =
      (  gCov[0][0]*gCov[1][1]*gCov[2][2]  
       + gCov[0][1]*gCov[1][2]*gCov[2][0]  
       + gCov[0][2]*gCov[1][0]*gCov[2][1]  
       - gCov[0][0]*gCov[1][2]*gCov[2][1]  
       - gCov[0][1]*gCov[1][0]*gCov[2][2]  
       - gCov[0][2]*gCov[1][1]*gCov[2][0])/gDet;
}

void geometry::setgCovInXCoords(const array XCoords[NDIM],
                                array gCov[NDIM][NDIM]
                               )
{
  switch (params::metric)
  {
    case metrics::MINKOWSKI:
      
      gCov[0][0] = -1.;
      gCov[1][1] = 1.;
      gCov[2][2] = 1.;
      gCov[3][3] = 1.;

      break;

    case metrics::MODIFIED_KERR_SCHILD:
      /* x^mu = {t, r, theta, phi}, X^mu = {t, X1, X2, phi} */

      array xCoords[NDIM];

      XCoordsToxCoords(XCoords, xCoords);

      /* Easier to read with (r, theta) than (x[1], x[2]) */
      array r     = xCoords[1];
      array theta = xCoords[2];

      /* r = exp(X1) => dr/dX = exp(X1) = r */
      array dr_dX1 = r;

      /* theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2) 
         => dtheta/dX2 = pi + pi*(1 - H_SLOPE)*cos(2*pi*X2) */
      array dtheta_dX2 = M_PI + M_PI*(1 - params::hSlope)*af::cos(2*M_PI*XCoords[2]);

      array sigma =  r*r + af::pow2(params::blackHoleSpin * af::cos(theta) );

      /* -(1 - 2*r/sigma) dt^2 */
      gCov[0][0] = -(1. - 2.*r/sigma);     

      /* (4*r/sigma * dr/dX1) dt dX1 */
      gCov[0][1] = (2.*r/sigma) * dr_dX1; 

      /* (0) dt dX2 */
      gCov[0][2] = 0.; 

      /* -(4*a*r*sin(theta)^2/sigma) dt dphi */
      gCov[0][3] = -(2.*params::blackHoleSpin*r*af::pow2(af::sin(theta))/sigma);

      /* (4*r/sigma * dr/dX1) dX1 dt */
      gCov[1][0] = gCov[0][1];

      /* ( (1 + 2*r/sigma)*dr/dX1*dr/dX1) dX1 dX1 */
      gCov[1][1] = (1. + 2*r/sigma) * dr_dX1 * dr_dX1;
  
      /* (0) dX1 dX2 */
      gCov[1][2] = 0.;

      /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
      gCov[1][3] =
        -params::blackHoleSpin*(1. + 2.*r/sigma)*af::pow2(af::sin(theta))*dr_dX1;

      /* (0) dX2 dt */
      gCov[2][0] = gCov[0][2];

      /* (0) dX2 dX1 */
      gCov[2][1] = gCov[1][2];

      /* (sigma*dtheta/dX2*dtheta/dX2) dX2 dX2 */
      gCov[2][2] = sigma*dtheta_dX2*dtheta_dX2;

      /* (0) dX2 dphi */
      gCov[2][3] = 0.;

      /* -(4*a*r*sin(theta)^2/sigma) dphi dt */
      gCov[3][0] = gCov[0][3];

      /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
      gCov[3][1] = gCov[1][3];

      /* (0) dphi dX2 */
      gCov[3][2] = gCov[2][3];

      /* (sin(theta)^2*(sigma + a^2*(1. + 2*r/sigma)*sin(theta)^2) dphi dphi */
      gCov[3][3] = (  af::pow2(af::sin(theta))
                    * (sigma +   af::pow2(params::blackHoleSpin*af::sin(theta))
                               * (1. + 2*r/sigma)
                      )
                   );
      
      break;
  }


}

void geometry::XCoordsToxCoords(const array XCoords[NDIM], 
                                array xCoords[NDIM]
                               )
{
  switch (params::metric)
  {
    case metrics::MINKOWSKI:
      
      for (int mu=0; mu<NDIM; mu++)
      {
        xCoords[mu] = XCoords[mu];
      }
      break;

    case metrics::MODIFIED_KERR_SCHILD:
      
      xCoords[0] = XCoords[0];
      xCoords[1] = af::exp(XCoords[1]);
      xCoords[2] =   M_PI*XCoords[2] 
                   + 0.5*(1 - params::hSlope)
                   * af::sin(2.*M_PI*XCoords[2]);
  
      xCoords[3] = XCoords[3];
      break;
  }
}


geometry::~geometry()
{
  delete xCoordsGrid;
  delete XCoordsGrid;
  delete alphaGrid;
  delete gGrid;
  delete gCovGrid;
  delete gConGrid;
  delete connectionGrid;
}


