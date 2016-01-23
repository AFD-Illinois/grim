#ifndef GRIM_GEOMETRY_H_
#define GRIM_GEOMETRY_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include <string>

void setXCoords(const grid &indices, int location, grid &XCoords);

class geometry
{
  private:
    array zero;
    array XCoords[3];
    void XCoordsToxCoords(const array XCoords[NDIM], array xCoords[NDIM]);
    void setgCovInXCoords(const array XCoords[NDIM], array gCov[NDIM][NDIM]);
    void setgDetAndgConFromgCov(const array gCov[NDIM][NDIM],
                                array &gDet, array gCon[NDIM][NDIM]
                               );
  public:
    int numGhost;

    array alpha;
    array g;
    array gCov[NDIM][NDIM];
    array gCon[NDIM][NDIM];

    array gammaUpDownDown[NDIM][NDIM][NDIM];

    geometry(const grid &XCoordsGrid);
    geometry(const geometry &geom,
             const int iStart, const int iEnd,
             const int jStart, const int jEnd,
             const int kStart, const int kEnd
            );
    ~geometry();


    void copyFrom(const geometry &geom,
                  const int iStart, const int iEnd,
                  const int jStart, const int jEnd,
                  const int kStart, const int kEnd
                 );
};

#endif /* GRIM_GEOMETRY_H_ */
