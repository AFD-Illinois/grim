#include "grim.h"
//#include <yaml-cpp/yaml.h>
#include <arrayfire.h>

inline int DELTA(int const &mu, int const &nu)
{
  return (mu==nu ? 1 : 0);
}

using af::array;
using af::span;
using af::where;

const int NDIM = 4;
const int LOCATIONS = 7;
const int AXISYM_LOCATIONS = 5;

namespace boundaries
{
  enum
  {
    PERIODIC, OUTFLOW, MIRROR, DIRICHLET
  };
};

namespace locations
{
  enum
  {
    CENTER, LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
  };
};

namespace directions
{
  enum
  {
    X1, X2, X3
  };
};

namespace timeStepping
{
  enum
  {
    EXPLICIT, IMEX, IMPLICIT
  };
};

namespace metrics
{
  enum
  {
    MINKOWSKI, MODIFIED_KERR_SCHILD
  };
};

namespace vars
{
  int RHO = 0;
  int U   = 1;
  int U1  = 2;
  int U2  = 3;
  int U3  = 4;
  int B1  = 5;
  int B2  = 6;
  int B3  = 7;
  int Q, DP;
  int dof;
};

namespace params
{
  int N1 = 8;
  int N2 = 8;
  int N3 = 8;
  int dim = 3;
  int numGhost = 2;

  int timeStepper = timeStepping::EXPLICIT;
  int metric = metrics::MINKOWSKI;
  double hSlope = 0.3;
  double blackHoleSpin = 0.9375;

  double X1Start = 0., X1End = 1.;
  double X2Start = 0., X2End = 1.;
  double X3Start = 0., X3End = 1.;

  int boundaryLeft   = boundaries::PERIODIC;
  int boundaryRight  = boundaries::PERIODIC;

  int boundaryTop    = boundaries::PERIODIC;
  int boundaryBottom = boundaries::PERIODIC;

  int boundaryFront  = boundaries::PERIODIC;
  int boundaryBack   = boundaries::PERIODIC;

  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;

  int conduction = 0;
  int viscosity  = 0;
  int highOrderTermsConduction = 1;
  int highOrderTermsViscosity = 1;
  double adiabaticIndex = 4./3;

  double slopeLimTheta = 2;

};

namespace gridParams
{
  bool haveGridParamsBeenSet = 0;
  int N1Local, N2Local, N3Local;
  int iLocalStart, jLocalStart, kLocalStart;
  int iLocalEnd,   jLocalEnd,   kLocalEnd;

  double dX1, dX2, dX3;
};

class grid
{
  public:
    DM dm, dmGhost;
    Vec globalVec, localVec;

    int numVars, numGhost;

    DMBoundaryType boundaryLeft,  boundaryRight;
    DMBoundaryType boundaryTop,   boundaryBottom;
    DMBoundaryType boundaryFront, boundaryBack;

    array *vars;

    grid(int numVars, int numGhost);
    ~grid();

    void setVarsWithVec(Vec vec);
};

class geometry
{
  private:
    array zero;
    void XCoordsToxCoords(const array XCoords[NDIM], array xCoords[NDIM]);
    void setgCovInXCoords(const array xCoords[NDIM], array gCov[NDIM][NDIM]);
    void setgDetAndgConFromgCov(const array gCov[NDIM][NDIM],
                                array gDet, array gCon[NDIM][NDIM]
                               );

  public:
    std::vector<int> axiSymmetricLocations
      = {locations::CENTER, 
         locations::LEFT, locations::RIGHT,
         locations::TOP,  locations::BOTTOM
        };

    std::vector<int> allLocations
      = {locations::CENTER, 
         locations::LEFT,  locations::RIGHT,
         locations::TOP,   locations::BOTTOM,
         locations::FRONT, locations::BACK
        };

    /* The grids manage allocations and IO */
    grid *xCoordsGrid;
    grid *XCoordsGrid;
    grid *alphaGrid, *gGrid;
    grid *gCovGrid, *gConGrid;
    grid *connectionGrid;

    /* Copy the data from the grids above into the arrays below so that they can
     * be conveniently used in physics routines */
    array xCoords[LOCATIONS][NDIM];
    array XCoords[LOCATIONS][NDIM];

    array alpha[AXISYM_LOCATIONS];
    array g[AXISYM_LOCATIONS];
    array gCov[AXISYM_LOCATIONS][NDIM][NDIM];
    array gCon[AXISYM_LOCATIONS][NDIM][NDIM];

    /* Connection coefficients only needed at CENTER */
    array gammaUpDownDown[NDIM][NDIM][NDIM];

    geometry(int numGhost);
    ~geometry();
};

class fluidElement
{
  private:
    const std::vector<std::vector<int>> indicesToLoopOver 
      = {{0},    {0, 1}, {0, 1}, {0, 2}, {0, 2}, {0, 3}, {0, 3}};
    /*   CENTER, LEFT,   RIGHT,  TOP,    BOTTOM, FRONT,  BACK*/

    array one;
  public:
    int loc;

    /* fluidElement parameters */
    array tau, chi, nu;
    
    array rho, u, u1, u2, u3, B1, B2, B3;
    array pressure, temperature;
    array qTilde, deltaPTilde;
    array q, deltaP;
  
    array gammaLorentzFactor, uCon[NDIM], uCov[NDIM];
    array bSqr, bCon[NDIM], bCov[NDIM];
    
    array NUp[NDIM];
    array TUpDown[NDIM][NDIM];

    fluidElement(const grid &prim, 
                 const geometry &geom, 
                 const int location
                );
    void set(const grid &prim, 
             const geometry &geom, 
             const int location
            );
    void setFluidElementParameters(const geometry &geom);
    void computeFluxes(const geometry &geom, 
                       const int direction,
                       grid &flux
                      );
    void computeSources(const geometry &geom,
                        grid &sources
                       );
};

geometry::geometry(int numGhost)
{
  xCoordsGrid     = new grid(LOCATIONS*NDIM,      numGhost);
  XCoordsGrid     = new grid(LOCATIONS*NDIM,      numGhost);

  alphaGrid       = new grid(AXISYM_LOCATIONS,           numGhost);
  gGrid           = new grid(AXISYM_LOCATIONS,           numGhost);
  gCovGrid        = new grid(AXISYM_LOCATIONS*NDIM*NDIM, numGhost);
  gConGrid        = new grid(AXISYM_LOCATIONS*NDIM*NDIM, numGhost);
  connectionGrid  = new grid(NDIM*NDIM*NDIM,             numGhost);

  grid indices(3, 0);

  /* Indices go from [-numGhost, N + numGhost) in each direction */
  indices.vars[directions::X1] = 
    af::range(af::dim4(gridParams::N1Local + 2*numGhost)) - numGhost;

  indices.vars[directions::X2] =
    af::range(af::dim4(gridParams::N2Local + 2*numGhost)) - numGhost;

  indices.vars[directions::X3] =
    af::range(af::dim4(gridParams::N3Local + 2*numGhost)) - numGhost;

  XCoordsGrid->vars[0 + NDIM*locations::CENTER] = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::RIGHT]  = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::LEFT]   = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::BOTTOM] = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::TOP]    = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::FRONT]  = 0.;
  XCoordsGrid->vars[0 + NDIM*locations::BACK]   = 0.;

  /* X1 coordinate at all the locations */
  XCoordsGrid->vars[1 + NDIM*locations::CENTER] = 
    params::X1Start + (indices.vars[directions::X1] + 0.5)*gridParams::dX1;

  XCoordsGrid->vars[1 + NDIM*locations::LEFT]   = 
    params::X1Start + (indices.vars[directions::X1]      )*gridParams::dX1;

  XCoordsGrid->vars[1 + NDIM*locations::RIGHT]  = 
    params::X1Start + (indices.vars[directions::X1] + 1. )*gridParams::dX1;

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
    params::X2Start + (indices.vars[directions::X2] + 0.5)*gridParams::dX2;

  XCoordsGrid->vars[2 + NDIM*locations::BOTTOM] = 
    params::X2Start + (indices.vars[directions::X2]      )*gridParams::dX2;

  XCoordsGrid->vars[2 + NDIM*locations::TOP]    = 
    params::X2Start + (indices.vars[directions::X2] + 1. )*gridParams::dX2;

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
    params::X3Start + (indices.vars[directions::X3] + 0.5)*gridParams::dX3;

  XCoordsGrid->vars[3 + NDIM*locations::BACK]   = 
    params::X3Start + (indices.vars[directions::X3]      )*gridParams::dX3;

  XCoordsGrid->vars[3 + NDIM*locations::FRONT]  = 
    params::X3Start + (indices.vars[directions::X3] + 1. )*gridParams::dX3;

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
  zero = af::constant(0., 
                      XCoords[locations::CENTER][0].dims(1),
                      XCoords[locations::CENTER][0].dims(2), 
                      XCoords[locations::CENTER][0].dims(3)
                     );


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
                                      array gDet,
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

class riemannSolver
{
  public:
    fluidElement *elemLeft, *elemRight;

    grid *primLeft, *primRight;
    grid *fluxLeft, *fluxRight;
    grid *consLeft, *consRight;

    riemannSolver(const geometry &geom);
    ~riemannSolver();

    void reconstruct(const grid &prim,
                     const int dir,
                     grid &primLeft,
                     grid &primRight
                    );

    void solve(const grid &prim,
               const geometry &geom,
               const int dir,
               grid &fluxes
              );
};

class timeStepper
{
  public:
    SNES snes;
    grid *prim, *primOld;
    grid *cons, *consOld;
    grid *sourcesOld;
    grid *residual;

    grid *fluxesX1, *fluxesX2, *fluxesX3;

    fluidElement *elem, *elemOld;
    riemannSolver *riemann;

    timeStepper(const geometry &geom);
    ~timeStepper();

    void timeStep();
};

PetscErrorCode computeResidual(SNES snes,
                               Vec primVec,
                               Vec residualVec,
                               void *ptr
                              )
{
  timeStepper *ts = (class timeStepper*)ptr;

}

timeStepper::timeStepper(const geometry &geom)
{

  prim        = new grid(vars::dof, params::numGhost);
  primOld     = new grid(vars::dof, params::numGhost);

  cons        = new grid(vars::dof, 0);
  consOld     = new grid(vars::dof, 0);

  sourcesOld  = new grid(vars::dof, 0);
  residual    = new grid(vars::dof, 0);

  switch (params::dim)
  {
    case 1:
      fluxesX1 = new grid(vars::dof, 1);
      break;

    case 2:
      fluxesX1 = new grid(vars::dof, 1);
      fluxesX2 = new grid(vars::dof, 1);
      break;

    case 3:
      fluxesX1 = new grid(vars::dof, 1);
      fluxesX2 = new grid(vars::dof, 1);
      fluxesX3 = new grid(vars::dof, 1);
      break;
  }

  elem    = new fluidElement(*prim, geom, locations::CENTER);
  elemOld = new fluidElement(*primOld, geom, locations::CENTER);

  riemann = new riemannSolver(geom);

  SNES snes;
  if (   params::timeStepper==timeStepping::EXPLICIT
      || params::timeStepper==timeStepping::IMEX
     )
  {
    SNESSetDM(snes, prim->dm);
  }
  else if (params::timeStepper==timeStepping::IMPLICIT)
  {
    SNESSetDM(snes, prim->dmGhost);
  }

  SNESSetFunction(snes, residual->globalVec, computeResidual, this);

}

timeStepper::~timeStepper()
{
  SNESDestroy(&snes);
  
  delete elem, elemOld;
  delete prim, primOld;
  delete cons, consOld;
  delete sourcesOld;
  delete residual;

  switch (params::dim)
  {
    case 1:
      delete fluxesX1;
      break;

    case 2:
      delete fluxesX1, fluxesX2;
      break;

    case 3:
      delete fluxesX1, fluxesX2, fluxesX3;
      break;
  }
}
//void timeStepper::timeStep()
//{
//  elemOld->set(primOld, geom, locations::CENTER);
//
//  elemOld->computeFluxes(geom, 0, consOld->vars);
//  elemOld->computeSources(geom, sourcesOld->vars);
////  
////  reconstruct(primOld->vars, 
////              primLeft->vars, primRight->vars,
////              directions::X1
////             );
////
//  riemann.solve(primRight->vars, primLeft->vars,
//                geom, directions::X1,
//                fluxesX1->vars
//               );
//
//}


int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscPrintf(PETSC_COMM_WORLD, 
              "#### System info ####\n"
             );
  if (rank==0)
  {
    af::info();
  };

  if      (params::conduction==1 && params::viscosity==0)
  {
    vars::Q   = 8;
    vars::dof = 9;
  }
  else if (params::conduction==0 && params::viscosity==1)
  {
    vars::DP  = 8;
    vars::dof = 9;
  }
  else if (params::conduction==1 && params::viscosity==1)
  {
    vars::Q   = 8;
    vars::DP  = 9;
    vars::dof = 10;
  }
  else if (params::conduction==0 && params::viscosity==0)
  {
    vars::dof = 8;
  }


  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    //timeStepper ts;
    
    grid prim(10, 0);
    grid flux(10, 0);
    geometry geom(0);
    for (int var=0; var<vars::dof; var++)
    {
      prim.vars[var] = af::randu(params::N1,
                                 params::N2,
                                 params::N3
                                );
    }


    fluidElement elem(prim, geom, locations::CENTER);
    elem.computeFluxes(geom, 0, flux);

    grid jacobian(10*10, 0);
    for (int column=0; column<vars::dof; column++)
    {
      primPlusEpsilon.vars[column] = prim.vars[column] + epsilon;
      
      elem.set(primPlusEpsilon, geom, locations::CENTER);

      for (int row=0; row<vars::dof; row++)
      {
        jacobian.vars[row + vars::dof*column]
          = conservedVars
      }
    }


  }

  PetscFinalize();  
  return(0);
}

grid::grid(int numVars, int numGhost)
{
  this->numVars = numVars;
  this->numGhost = numGhost;

  /* Implementations for MIRROR, OUTFLOW in boundary.cpp and DIRICHLET in
   * problem.cpp */
  boundaryLeft  = DM_BOUNDARY_GHOSTED; boundaryRight  = DM_BOUNDARY_GHOSTED;
  boundaryTop   = DM_BOUNDARY_GHOSTED; boundaryBottom = DM_BOUNDARY_GHOSTED;
  boundaryFront = DM_BOUNDARY_GHOSTED; boundaryBack   = DM_BOUNDARY_GHOSTED;

  if (   params::boundaryLeft  == boundaries::PERIODIC 
      || params::boundaryRight == boundaries::PERIODIC
     )
  {
    boundaryLeft  = DM_BOUNDARY_PERIODIC;
    boundaryRight = DM_BOUNDARY_PERIODIC;
  }

  if (   params::boundaryTop    == boundaries::PERIODIC 
      || params::boundaryBottom == boundaries::PERIODIC
     )
  {
    boundaryTop    = DM_BOUNDARY_PERIODIC;
    boundaryBottom = DM_BOUNDARY_PERIODIC;
  }

  if (   params::boundaryFront == boundaries::PERIODIC 
      || params::boundaryBack  == boundaries::PERIODIC
     )
  {
    boundaryFront = DM_BOUNDARY_PERIODIC;
    boundaryBack  = DM_BOUNDARY_PERIODIC;
  }

  switch (params::dim)
  {
    case 1:
      DMDACreate1d(PETSC_COMM_WORLD, boundaryLeft, 
                   params::N1, numVars, 0, NULL,
                   &dm
                  );

      DMDACreate1d(PETSC_COMM_WORLD, boundaryLeft, 
                   params::N1, numVars, numGhost, NULL,
                   &dmGhost
                  );

      break;
  
    case 2:
      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, 0,
                   PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      DMDACreate2d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost,
                   PETSC_NULL, PETSC_NULL,
                   &dmGhost
                  );
      break;

    case 3:
      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2, params::N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, 0, 
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      DMDACreate3d(PETSC_COMM_WORLD, 
                   boundaryLeft, boundaryBottom, boundaryBack,
                   DMDA_STENCIL_BOX,
                   params::N1, params::N2, params::N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhost, 
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dmGhost
                  );
      break;
  }

  DMCreateGlobalVector(dm, &globalVec);
  DMCreateLocalVector(dmGhost, &localVec);

  if (!gridParams::haveGridParamsBeenSet)
  {
    DMDAGetCorners
      (dm,
       &gridParams::iLocalStart,
       &gridParams::jLocalStart,
       &gridParams::kLocalStart,
       &gridParams::N1Local,
       &gridParams::N2Local,
       &gridParams::N3Local
      );

    gridParams::iLocalEnd = gridParams::iLocalStart + gridParams::N1Local;
    gridParams::jLocalEnd = gridParams::jLocalStart + gridParams::N2Local;
    gridParams::kLocalEnd = gridParams::kLocalStart + gridParams::N3Local;

    gridParams::dX1 = (params::X1End - params::X1Start)/params::N1;
    gridParams::dX2 = (params::X2End - params::X2Start)/params::N2;
    gridParams::dX3 = (params::X3End - params::X3Start)/params::N3;
  }

  array varsCopiedFromVec;
  if (numGhost > 0)
  {
    DMGlobalToLocalBegin(dmGhost, globalVec, INSERT_VALUES, localVec);
    DMGlobalToLocalEnd(dmGhost, globalVec, INSERT_VALUES, localVec);

    double *pointerToLocalVec;
    VecGetArray(localVec, &pointerToLocalVec);

    switch (params::dim)
    {
      case 1:
        varsCopiedFromVec = 
          array(numVars,
                gridParams::N1Local + 2*numGhost, 
                gridParams::N2Local,
                gridParams::N3Local,
                pointerToLocalVec
               );
        break;

      case 2:
        varsCopiedFromVec = 
          array(numVars,
                gridParams::N1Local + 2*numGhost, 
                gridParams::N2Local + 2*numGhost,
                gridParams::N3Local,
                pointerToLocalVec
               );
        break;

      case 3:
        varsCopiedFromVec = 
          array(numVars,
                gridParams::N1Local + 2*numGhost, 
                gridParams::N2Local + 2*numGhost,
                gridParams::N3Local + 2*numGhost,
                pointerToLocalVec
               );
        break;
      }

    VecRestoreArray(localVec, &pointerToLocalVec);
  }
  else
  {
    double *pointerToGlobalVec;
    VecGetArray(globalVec, &pointerToGlobalVec);

    varsCopiedFromVec = 
      array(numVars,
            gridParams::N1Local, gridParams::N2Local, gridParams::N3Local,
            pointerToGlobalVec
           );

    VecRestoreArray(globalVec, &pointerToGlobalVec);
  }

  vars = new array[numVars];
  for (int var=0; var<numVars; var++)
  {
    vars[var] = varsCopiedFromVec(var, span, span, span);
  }

}

void grid::setVarsWithVec(Vec vec)
{
  array varsCopiedFromVec;
  if (numGhost > 0)
  {
    DMGlobalToLocalBegin(dmGhost, vec, INSERT_VALUES, localVec);
    DMGlobalToLocalEnd(dmGhost, vec, INSERT_VALUES, localVec);

    double *pointerToLocalVec;
    VecGetArray(localVec, &pointerToLocalVec);

    varsCopiedFromVec = 
      array(numVars,
            gridParams::N1Local + 2*numGhost, 
            gridParams::N2Local + 2*numGhost,
            gridParams::N3Local + 2*numGhost,
            pointerToLocalVec
           );

    VecRestoreArray(globalVec, &pointerToLocalVec);
  }
  else
  {
    double *pointerToGlobalVec;
    VecGetArray(vec, &pointerToGlobalVec);

    varsCopiedFromVec = 
      array(numVars,
            gridParams::N1Local, gridParams::N2Local, gridParams::N3Local,
            pointerToGlobalVec
           );

    VecRestoreArray(vec, &pointerToGlobalVec);
  }

  for (int var=0; var<vars::dof; var++)
  {
    vars[var] = varsCopiedFromVec(var, span, span, span);
  }
}

grid::~grid()
{
  delete [] vars;
  VecDestroy(&globalVec);
  VecDestroy(&localVec);

  DMDestroy(&dm);
  DMDestroy(&dmGhost);
}

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

riemannSolver::riemannSolver(const geometry &geom)
{
  primLeft  = new grid(vars::dof, params::numGhost);
  primRight = new grid(vars::dof, params::numGhost);

  elemLeft  = new fluidElement(*primLeft, geom, locations::LEFT);
  elemRight = new fluidElement(*primRight, geom, locations::LEFT);

  fluxLeft  = new grid(vars::dof, 0);
  fluxRight = new grid(vars::dof, 0);

  consLeft  = new grid(vars::dof, 0);
  consRight = new grid(vars::dof, 0);
}

riemannSolver::~riemannSolver()
{
  delete primLeft, primRight;
  delete elemLeft, elemRight;
  delete fluxLeft, fluxRight;
  delete consLeft, consRight;
}

void riemannSolver::solve(const grid &prim,
                          const geometry &geom,
                          const int dir,
                          grid &flux
                         )
{
  int location;
  switch (dir)
  {
    case directions::X1:
      location = locations::LEFT;
      break;

    case directions::X2:
      location = locations::BOTTOM;
      break;

    case directions::X3:
      location = locations::BACK;
      break;
  }

  reconstruct(prim, dir, *primLeft, *primRight);

  elemLeft->set(*primRight, geom, location);
  elemRight->set(*primLeft, geom, location);

  elemLeft->computeFluxes(geom, dir, *fluxLeft);
  elemLeft->computeFluxes(geom, 0,   *consLeft);

  elemRight->computeFluxes(geom, dir, *fluxRight);
  elemRight->computeFluxes(geom, 0,   *consRight);

  double cLaxFriedrichs = 1.;

  for (int var=0; var<vars::dof; var++)
  {
    flux.vars[var] = 0.5*(  fluxLeft->vars[var] + fluxRight->vars[var]
                          - cLaxFriedrichs 
                          * (consRight->vars[var] - consLeft->vars[var])
                         );
  }
            
}

array minmod(array x, array y, array z)
{
  array minOfAll = af::min(
                           af::min(af::abs(x), af::abs(y)), 
                           af::abs(z)
                          );

  return 0.25 * abs(sign(x) + sign(y) ) * (sign(x) + sign(z) ) * minOfAll;
}

void riemannSolver::reconstruct(const grid &prim,
                                const int dir,
                                grid &primLeft,
                                grid &primRight
                               )
{
  float filter1D[]  = {1,-1, 0, /* Forward difference */
                       0, 1,-1  /* Backward difference */
                      };
  array filter;
  switch (dir)
  {
    case directions::X1:
      filter =  array(3, 1, 1, 2, filter1D);
      break;

    case directions::X2:
      filter =  array(1, 3, 1, 2, filter1D);
      break;

    case directions::X3:
      filter =  array(1, 1, 3, 2, filter1D);
      break;
  }

  for(int var=0; var<vars::dof; var++)
  {
    array dvar_dX = convolve(prim.vars[var], filter);

    array forwardDiff  = dvar_dX(span, span, span, 0);
    array backwardDiff = dvar_dX(span, span, span, 1);
    array centralDiff  = backwardDiff + forwardDiff;

    array slope =  minmod(params::slopeLimTheta * backwardDiff, 
                          0.5 * centralDiff, 
                          params::slopeLimTheta * forwardDiff
                         );

    primLeft.vars[var]  = prim.vars[var] - 0.5*slope;
    primRight.vars[var] = prim.vars[var] + 0.5*slope;
  }

}


                                 
