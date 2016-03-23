from gridHeaders cimport grid

cdef extern from "boundary.hpp":
  cdef enum:
    BOUNDARIES_PERIODIC  "boundaries::PERIODIC"
    BOUNDARIES_MIRROR    "boundaries::MIRROR"
    BOUNDARIES_OUTFLOW   "boundaries::OUTFLOW"
    BOUNDARIES_DIRICHLET "boundaries::DIRICHLET"

cdef extern from "boundary.hpp" namespace "boundaries":
  void applyBoundaryConditions(const int boundaryLeft, const int boundaryRight,
                               const int boundaryTop, const int boundaryBottom,
                               const int boundaryFront, const int boundaryBack,
                               grid &prim
                              )
