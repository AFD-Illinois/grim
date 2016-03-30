"""
timeStepperPy
======

Python interface to the timeStepper class.
"""

import numpy as np
from gridPy cimport gridPy, coordinatesGridPy
from geometryPy cimport geometryPy
from physicsPy cimport fluidElementPy
from timestepperHeaders cimport timeStepper
from timestepperHeaders cimport STAGE_HALF_STEP
from timestepperHeaders cimport STAGE_FULL_STEP

# Time stepping stage macros
HALF_STEP = STAGE_HALF_STEP
FULL_STEP = STAGE_FULL_STEP

cdef class timeStepperPy(object):

  def __cinit__(self, const int N1,
                      const int N2,
                      const int N3,
                      const int dim,
                      const int numVars,
                      const int numGhost,
                      const double time,
                      const double dt,
                      const int boundaryLeft, const int boundaryRight,
                      const int boundaryTop,  const int boundaryBottom,
                      const int boundaryFront, const int boundaryBack,
                      const int metric,
                      const double blackHoleSpin,
                      const double hSlope,
                      const double X1Start, const double X1End,
                      const double X2Start, const double X2End,
                      const double X3Start, const double X3End
               ):
    self.timeStepperPtr = \
        new timeStepper(N1, N2, N3, dim, numVars, numGhost,
                        time, dt,
                        boundaryLeft, boundaryRight,
                        boundaryTop,  boundaryBottom,
                        boundaryFront, boundaryBack,
                        metric, blackHoleSpin, hSlope,
                        X1Start, X1End,
                        X2Start, X2End,
                        X3Start, X3End
                       )

    self.XCoords = coordinatesGridPy()
    self.XCoords.setGridPtr(self.timeStepperPtr.XCoords)

    self.prim         = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.prim)
    self.primOld      = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.primOld)
    self.primHalfStep = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.primHalfStep)

    self.fluxesX1 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.fluxesX1)
    self.fluxesX2 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.fluxesX2)
    self.fluxesX3 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.fluxesX3)

    self.sourcesExplicit = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.sourcesExplicit)

    self.divB = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.divB)

    self.geomCenter = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomCenter)
    self.geomLeft = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomLeft)
    self.geomRight = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomRight)
    self.geomBottom = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomBottom)
    self.geomTop = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomTop)

    self.elem = \
        fluidElementPy.createFluidElementPyFromElemPtr(self.timeStepperPtr.elem)


  def __dealloc__(self):
    del self.timeStepperPtr

  property XCoords:
    def __get__(self):
     return self.XCoords

  property elem:
    def __get__(self):
     return self.elem

  property geomCenter:
    def __get__(self):
     return self.geomCenter

  property geomRight:
    def __get__(self):
     return self.geomRight

  property geomLeft:
    def __get__(self):
     return self.geomLeft

  property geomBottom:
    def __get__(self):
     return self.geomBottom

  property geomTop:
    def __get__(self):
     return self.geomTop

  property prim:
    def __get__(self):
     return self.prim

  property primOld:
    def __get__(self):
     return self.primOld

  property primHalfStep:
    def __get__(self):
     return self.primHalfStep

  property fluxesX1:
    def __get__(self):
     return self.fluxesX1

  property fluxesX2:
    def __get__(self):
     return self.fluxesX2

  property fluxesX3:
    def __get__(self):
     return self.fluxesX3

  property sourcesExplicit:
    def __get__(self):
     return self.sourcesExplicit

  property divB:
    def __get__(self):
     return self.divB


  def fluxCT(self):
     cdef int numReads  = 0
     cdef int numWrites = 0
     self.timeStepperPtr.fluxCT(numReads, numWrites)
     return numReads, numWrites

  def computeDivB(self, gridPy prim):
     cdef int numReads  = 0
     cdef int numWrites = 0
     self.timeStepperPtr.computeDivB(prim.getGridPtr()[0],
                                     numReads, numWrites
                                    )
     return numReads, numWrites

  def timeStep(self):
      cdef int numReads  = 0
      cdef int numWrites = 0
      self.timeStepperPtr.timeStep(numReads, numWrites)
      return numReads, numWrites

  def computeDivOfFluxes(self, gridPy prim):
    cdef int numReads  = 0
    cdef int numWrites = 0
    self.timeStepperPtr.computeDivOfFluxes(prim.getGridPtr()[0],
                                           numReads, numWrites
                                          )
    return numReads, numWrites
