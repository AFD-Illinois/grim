#include "radiation.hpp"

emission::emission(const fluidElement &elem_, photons &ph_):
  ph(ph_), /* Initialization list */
  elem(elem_)
{
}

emission::~emission()
{
  /* Nothing to do here */
}

void emission::emit(const double dt)
{
}
