#ifndef GRIM_RADIATION_H_
#define GRIM_RADIATION_H_

#include "../../params.hpp"
#include "../../grid/grid.hpp"
#include "../../geometry/geometry.hpp"
#include <algorithm> // std::remove_if

class photonManager
{
  public:
    int numPhotons;
    //int photonSize;

    /* Radiation four-force vector */

    photonManager(int numPhotons);
    ~photonManager();

    void addPhotons(const int numToAdd);
    void updateOne();
    void update(const double dt);

    std::vector<af::array> X;
    std::vector<af::array> K;
    std::vector<af::array> Xhalf;
    std::vector<af::array> Khalf;
    std::vector<af::array> dK;
    
    af::array dl;
    af::array w;
    af::array E;
    af::array K0Init;
    af::array K3Init;
    af::array nScatt;
    af::array iOrigin;
    af::array jOrigin;
    af::array isActive;

};

#endif /* GRIM_RADIATION_H_ */
