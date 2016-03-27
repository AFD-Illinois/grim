"""
timeStepperPy
======

Python interface to the timeStepper class.
"""

import numpy as np
from gridPy cimport gridPy, coordinatesGridPy
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


  def __dealloc__(self):
    del self.timeStepperPtr

  property XCoords:
    def __get__(self):
     return self.XCoords

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

  def fluxCT(self):
     cdef int numReads  = 0
     cdef int numWrites = 0
     self.timeStepperPtr.fluxCT(numReads, numWrites)
     return numReads, numWrites

  def computeDivB(self, gridPy prim, gridPy divB):
     cdef int numReads  = 0
     cdef int numWrites = 0
     self.timeStepperPtr.computeDivB(prim.getGridPtr()[0],
                                     divB.getGridPtr()[0],
                                     numReads, numWrites
                                    )
     return numReads, numWrites

  def timeStep(self):
      cdef int numReads  = 0
      cdef int numWrites = 0
      self.timeStepperPtr.timeStep(numReads, numWrites)
      return numReads, numWrites
