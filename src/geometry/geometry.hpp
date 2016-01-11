#ifndef GRIM_GEOMETRY_H_
#define GRIM_GEOMETRY_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include <string>

class geometry
{
  private:
    array zero;
    void XCoordsToxCoords(const array XCoords[NDIM], array xCoords[NDIM]);
    void setgCovInXCoords(const array xCoords[NDIM], array gCov[NDIM][NDIM]);
    void setgDetAndgConFromgCov(const array gCov[NDIM][NDIM],
                                array &gDet, array gCon[NDIM][NDIM]
                               );

  public:
    int numGhost;

    std::vector<int> allLocations
      = {locations::CENTER, 
         locations::LEFT,  locations::RIGHT,
         locations::TOP,   locations::BOTTOM,
         locations::FRONT, locations::BACK
        };

    grid *XCoordsGrid;

    array xCoords[LOCATIONS][NDIM];
    array XCoords[LOCATIONS][NDIM];

    array alpha[LOCATIONS];
    array g[LOCATIONS];
    array gCov[LOCATIONS][NDIM][NDIM];
    array gCon[LOCATIONS][NDIM][NDIM];

    /* Connection coefficients only needed at CENTER */
    array gammaUpDownDown[NDIM][NDIM][NDIM];

    geometry();
    ~geometry();

//    void getHostPtrTo(std::string str);
};


#endif /* GRIM_GEOMETRY_H_ */
