#ifndef GRIM_GEOMETRY_H_
#define GRIM_GEOMETRY_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include <string>

void setXCoords(const int location, grid &XCoords);

class geometry
{
  private:
    array zero;
    array XCoords[3];
    void setgCovInXCoords(const array XCoords[NDIM], array gCov[NDIM][NDIM]);
    void setgDetAndgConFromgCov(const array gCov[NDIM][NDIM],
                                array &gDet, array gCon[NDIM][NDIM]
                               );
    void computeGammaDownDownDown(const int eta,
                                  const int mu,
                                  const int nu,
                                  array& out
                                 );
  public:
    int numGhost;

    array alpha;
    array g;
    array gCov[NDIM][NDIM];
    array gCon[NDIM][NDIM];

    array gammaUpDownDown[NDIM][NDIM][NDIM];

    geometry(const grid &XCoordsGrid);
    ~geometry();

    void computeConnectionCoeffs();
    void getXCoords(array tXCoords[3]) const
    {
      for(int d=0;d<3;d++) tXCoords[d]=XCoords[d];
    }
    void XCoordsToxCoords(const array XCoords[3], array xCoords[3]) const;
};

#endif /* GRIM_GEOMETRY_H_ */
