#include "geometry.hpp"

//void geometry::getHostPtrTo(std::string str)
//{
//  if (str=="xCoords")
//  {
//    xCoordsHostPtr = 
//
//
//  }
//
//}

geometry::geometry()
{
  numGhost  = params::numGhost;

  XCoordsGrid = new grid(LOCATIONS*NDIM);

  /* Indices go from [-numGhost, N + numGhost) in each direction */
  array indicesX1 = 
    af::range(XCoordsGrid->N1Total, /* number of total zones in X1 */
              XCoordsGrid->N2Total, /* number of total zones in X2 */
              XCoordsGrid->N3Total, /* number of total zones in X3 */
              1,                    /* number of variables */
              directions::X1 ,      /* Vary indices in X1 direction */
              f64                   /* Double precision */
             ) - numGhost;

  array indicesX2 = 
    af::range(XCoordsGrid->N1Total, /* number of total zones in X1 */
              XCoordsGrid->N2Total, /* number of total zones in X2 */
              XCoordsGrid->N3Total, /* number of total zones in X3 */
              1,                    /* number of variables */
              directions::X2 ,      /* Vary indices in X2 direction */
              f64                   /* Double precision */
             ) - numGhost;

  array indicesX3 = 
    af::range(XCoordsGrid->N1Total, /* number of total zones in X1 */
              XCoordsGrid->N2Total, /* number of total zones in X2 */
              XCoordsGrid->N3Total, /* number of total zones in X3 */
              1,                    /* number of variables */
              directions::X3 ,      /* Vary indices in X3 direction */
              f64                   /* Double precision */
             ) - numGhost;

  /* Use the following array to allocate other arrays */
  zero = 0.*XCoordsGrid->vars[0];

  /* Temporal X coordinate. Just set to zero */
  for (int loc: allLocations)
  {
    XCoords[loc][0] = zero;
  }

  double dX1 = XCoordsGrid->dX1;
  double dX2 = XCoordsGrid->dX2;
  double dX3 = XCoordsGrid->dX3;

  /* X1 coordinate at all the locations */
  XCoords[locations::CENTER][1] = params::X1Start + (indicesX1 + 0.5)*dX1;
  XCoords[locations::LEFT][1]   = params::X1Start + (indicesX1      )*dX1;
  XCoords[locations::RIGHT][1]  = params::X1Start + (indicesX1 + 1. )*dX1;
  XCoords[locations::BOTTOM][1] = XCoords[locations::CENTER][1];
  XCoords[locations::TOP][1]    = XCoords[locations::CENTER][1];
  XCoords[locations::FRONT][1]  = XCoords[locations::CENTER][1];
  XCoords[locations::BACK][1]   = XCoords[locations::CENTER][1];

  /* X2 coordinate at all the locations */
  XCoords[locations::CENTER][2] = params::X2Start + (indicesX2 + 0.5)*dX2;
  XCoords[locations::LEFT][2]   = XCoords[locations::CENTER][2];
  XCoords[locations::RIGHT][2]  = XCoords[locations::CENTER][2];
  XCoords[locations::BOTTOM][2] = params::X2Start + (indicesX2      )*dX2;
  XCoords[locations::TOP][2]    = params::X2Start + (indicesX2 + 1. )*dX2;
  XCoords[locations::FRONT][2]  = XCoords[locations::CENTER][2];
  XCoords[locations::BACK][2]   = XCoords[locations::CENTER][2];

  /* X3 coordinate at all the locations */
  XCoords[locations::CENTER][3] = params::X3Start + (indicesX3 + 0.5)*dX3;
  XCoords[locations::LEFT][3]   = XCoords[locations::CENTER][3];
  XCoords[locations::RIGHT][3]  = XCoords[locations::CENTER][3];
  XCoords[locations::BOTTOM][3] = XCoords[locations::CENTER][3];
  XCoords[locations::TOP][3]    = XCoords[locations::CENTER][3];
  XCoords[locations::FRONT][3]  = params::X3Start + (indicesX3 + 1. )*dX3;
  XCoords[locations::BACK][3]   = params::X3Start + (indicesX3      )*dX3;

  /* Compute xCoords from XCoords */
  for (int loc: allLocations)
  {
    XCoordsToxCoords(XCoords[loc], xCoords[loc]);
  }

  for (int loc : allLocations)
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
  }

  /* Only care about connection coefficients at the center */
  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int lam = 0; lam<NDIM; lam++)
      {
      	gammaUpDownDown[mu][nu][lam]=zero;
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
  delete XCoordsGrid;
}


