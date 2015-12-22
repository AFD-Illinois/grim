#ifndef GRIM_GEOMETRY_H_
#define GRIM_GEOMETRY_H_

#include "../params.hpp"
#include "../grid/grid.hpp"

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


#endif /* GRIM_GEOMETRY_H_ */
