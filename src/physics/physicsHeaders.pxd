from gridHeaders cimport grid, coordinatesGrid
from geometryHeaders cimport geometry

cdef extern from "physics.hpp":
  cdef cppclass fluidElement:
    fluidElement(const grid &prim,
                 const geometry &geom,
                 int &numReads,
                 int &numWrites
                )
    void computeFluxes(const int direction,
                       grid &flux,
                       int &numReads,
                       int &numWrites
                      )
