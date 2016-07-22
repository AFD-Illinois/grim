#ifndef GRIM_RADIATION_H_
#define GRIM_RADIATION_H_

#include "../../params.hpp"
#include "../../grid/grid.hpp"
#include "../../geometry/geometry.hpp"
#include "../physics.hpp"
#include <algorithm> // std::remove_if

/* The photons data structure. Handles what the photons can do
 * themselves, namely move along their geodesics */
class photons
{
  public:
    int numPhotons;

    /* Radiation four-force vector */

    photons(const geometry& geom_, const int numPhotons);
    ~photons();

    /* This functionality provided by emission and interactions */
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
    
  private:
    const geometry& geom;

};

/* The home of emission processes. */
class emission
{
  public:
    emission(const fluidElement &elem_, photons &ph_);
    ~emission();
    
    void emit(const double dt);
    
  private:
    photons& ph;
    const fluidElement& elem;
    
};

#endif /* GRIM_RADIATION_H_ */
