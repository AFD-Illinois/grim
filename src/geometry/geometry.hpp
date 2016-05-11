#ifndef GRIM_GEOMETRY_H_
#define GRIM_GEOMETRY_H_

#include "../params.hpp"
#include "../grid/grid.hpp"

class geometry
{
  private:
    array zero;
    array XCoords[3];
    array xCoords[3];
    void setgCovInXCoords(const array XCoords[NDIM], array gCov[NDIM][NDIM]);
    void setgDetAndgConFromgCov(const array gCov[NDIM][NDIM],
                                array &gDet, array gCon[NDIM][NDIM]
                               );
    void computeGammaDownDownDown(const int eta,
                                  const int mu,
                                  const int nu,
                                  array& out,
				  const array gCovPlus[NDIM-1][NDIM][NDIM],
				  const array gCovMinus[NDIM-1][NDIM][NDIM]
                                 );
    // Change in coord value when computing metric derivatives
    double GAMMA_EPS;
  
  public:
    int N1, N2, N3, dim, numGhost;

    int metric;
    double blackHoleSpin;
    double hSlope;

    array alpha;
    array g;
    array gCov[NDIM][NDIM];
    array gCon[NDIM][NDIM];

    array gammaUpDownDown[NDIM][NDIM][NDIM];

    geometry(const int metric,
             const double blackHoleSpin,
             const double hSlope,
             const coordinatesGrid &XCoordsGrid
            );
    ~geometry();

    void computeConnectionCoeffs();
    void getXCoords(array tXCoords[3]) const
    {
      for(int d=0;d<3;d++) tXCoords[d]=XCoords[d];
    }
    void getxCoords(array txCoords[3]) const
    {
      for(int d=0;d<3;d++) txCoords[d]=xCoords[d];
    }
    void XCoordsToxCoords(const array XCoords[3], array xCoords[3]) const;

    /* Pointers to data on host. Needed to get data into Numpy */
    grid *gCovGrid;
    grid *gConGrid;
    grid *gGrid;
    grid *alphaGrid;
    grid *gammaUpDownDownGrid;
    grid *xCoordsGrid;

    void setgConGrid();
    void setgCovGrid();
    void setgGrid();
    void setalphaGrid();
    void setgammaUpDownDownGrid();
    void setxCoordsGrid();
};

#endif /* GRIM_GEOMETRY_H_ */
