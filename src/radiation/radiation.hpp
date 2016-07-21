#ifndef GRIM_RADIATION_H_
#define GRIM_RADIATION_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../geometry/geometry.hpp"

class photonManager
{
  public:
    std::vector photons;

    photonManager();
    ~photonManager();

    void addPhoton(photon *ph);

};

#endif /* GRIM_RADIATION_H_ */
