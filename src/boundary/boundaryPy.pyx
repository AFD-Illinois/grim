"""
boundaryPy
======

Python interface to the boundary conditions
"""
from gridPy cimport gridPy

from boundaryHeaders cimport BOUNDARIES_OUTFLOW
from boundaryHeaders cimport BOUNDARIES_MIRROR
from boundaryHeaders cimport BOUNDARIES_DIRICHLET
from boundaryHeaders cimport BOUNDARIES_PERIODIC

from boundaryHeaders cimport applyBoundaryConditions

# boundary macros
OUTFLOW   = BOUNDARIES_OUTFLOW
MIRROR    = BOUNDARIES_MIRROR
DIRICHLET = BOUNDARIES_DIRICHLET
PERIODIC  = BOUNDARIES_PERIODIC

def applyBoundaryConditionsPy(int boundaryLeft,  int boundaryRight,
                              int boundaryTop,   int boundaryBottom,
                              int boundaryFront, int boundaryBack,
                              gridPy prim
                             ):
  applyBoundaryConditions(boundaryLeft,  boundaryRight,
                          boundaryTop,   boundaryBottom,
                          boundaryFront, boundaryBack,
                          prim.getGridPtr()[0]
                         )
