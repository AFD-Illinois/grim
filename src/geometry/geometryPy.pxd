import numpy as np
cimport numpy as np
from geometryHeaders cimport geometry

cdef class geometryPy(object):
  cdef geometry *geometryPtr
  cdef np.ndarray gCov
  cdef np.ndarray gCon
  cdef np.ndarray g
  cdef np.ndarray alpha
  cdef np.ndarray gammaUpDownDown
  cdef np.ndarray xCoords
