import numpy as np
cimport numpy as np
from timestepperHeaders cimport timeStepper
from gridPy cimport gridPy, coordinatesGridPy

cdef class timeStepperPy(object):
  cdef timeStepper *timeStepperPtr
  cdef coordinatesGridPy XCoords
  cdef gridPy prim, primOld, primHalfStep
  cdef gridPy fluxesX1, fluxesX2, fluxesX3
