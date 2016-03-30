import numpy as np
cimport numpy as np
from timestepperHeaders cimport timeStepper
from gridPy cimport gridPy, coordinatesGridPy
from geometryPy cimport geometryPy
from physicsPy cimport fluidElementPy

cdef class timeStepperPy(object):
  cdef timeStepper *timeStepperPtr
  cdef coordinatesGridPy XCoords
  cdef geometryPy geomCenter
  cdef geometryPy geomLeft, geomRight
  cdef geometryPy geomTop, geomBottom
  cdef geometryPy geomFront, geomBack
  cdef gridPy prim, primOld, primHalfStep
  cdef gridPy fluxesX1, fluxesX2, fluxesX3
  cdef gridPy sourcesExplicit
  cdef gridPy divB
  cdef fluidElementPy elem, elemOld, elemHalfStep
