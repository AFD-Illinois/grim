#include "radiation.hpp"

photons::photons(const geometry& geom_, const int numPhotons_):
  geom(geom_)
{
  numPhotons = numPhotons_;
  
  /* Create pools of memory for photons */
  for (int mu = 0; mu < NDIM; mu++)
  {
    X.push_back(    af::constant(0, numPhotons, f64));
    K.push_back(    af::constant(0, numPhotons, f64));
    Xhalf.push_back(af::constant(0, numPhotons, f64));
    Khalf.push_back(af::constant(0, numPhotons, f64));
    dK.push_back(   af::constant(0, numPhotons, f64));
  }
    
  dl       = af::constant(0, numPhotons, f64);
  w        = af::constant(0, numPhotons, f64);
  E        = af::constant(0, numPhotons, f64);
  K0Init   = af::constant(0, numPhotons, f64);
  K3Init   = af::constant(0, numPhotons, f64);
  nScatt   = af::constant(0, numPhotons, u32);
  iOrigin  = af::constant(0, numPhotons, u32);
  jOrigin  = af::constant(0, numPhotons, u32);
  isActive = af::constant(0, numPhotons, u16);
}

photons::~photons()
{
  /* Nothing to do here? */
}

/* Add an array of photons to manager's pool */
void photons::addPhotons(const int numToAdd)
{
  /* Ensure that pool does not overflow */
  if (numToAdd + af::where(isActive > 0).elements() > numPhotons)
  {
    PetscPrintf(PETSC_COMM_WORLD, "Trying to create more photons than pool allows!\n");
    PetscPrintf(PETSC_COMM_WORLD, "Tried to create:  %d photons\n", 
                numToAdd);
    PetscPrintf(PETSC_COMM_WORLD, "Currently active: %d photons\n", 
                af::where(isActive > 0).elements());
    PetscPrintf(PETSC_COMM_WORLD, "Pool size:        %d photons\n", 
                numPhotons);
    exit(-1);
  } 
  else /* Add photons to manager pool */ 
  {
    /* Get indices of inactive elements in pool */
    array inactive = af::where(isActive == 0);
    array needed = inactive(af::seq(0, numToAdd-1));
    isActive(needed) = 1;
  }
  
}

/* Update a photon */
void photons::updateOne()
{
  
}

/*bool erasePredicate()
{
  return 0;
}*/

/* Update all photons by dt */
void photons::update(const double dt)
{
  /* Integrate geodesics of active photons by dt */
  int numActive = af::where(isActive > 0).elements();
  if(numActive > 0)
  {
    printf("%d living particles\n", numActive);
    for (int mu = 0; mu < NDIM; mu++)
    {
      // X[mu](isActive > 0) += ... 
    }
  }
  
  /* Remove any dead photons from memory */
  /*photons.erase(std::remove_if(photons.begin(), 
                               photons.end(),
                               erasePredicate),
                               photons.end());*/
}

