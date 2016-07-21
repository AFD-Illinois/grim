#include "radiation.hpp"

void photonManager::photonManager()
{
}

/* Add a photon to manager's ledger */
void photonManager::addPhoton(photon *ph)
{
  photons.push_back(ph);
}

/* Update photons and record whether they should be removed */
bool photonManager::updateOne(const double dt, photon *ph)
{
  return ph->push(dt);
}

/* Update all photons by dt */
void photonManager::update(const double dt)
{
  photons.erase(std::remove_if(photons.begin(), 
                               photons.end(),
                               updateOne),
                               photons.end());
}

