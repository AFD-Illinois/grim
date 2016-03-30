import numpy as np
cimport numpy as np
from geometryHeaders cimport geometry

cdef class geometryPy(object):
  cdef geometry *geometryPtr
  cdef int usingExternalPtr
  cdef np.ndarray gCov
  cdef np.ndarray gCon
  cdef np.ndarray g
  cdef np.ndarray alpha
  cdef np.ndarray gammaUpDownDown
  cdef np.ndarray xCoords
  cdef geometry* getGeometryPtr(self)
  cdef setGeometryPtr(self, geometry *geometryPtr)

  @staticmethod
  cdef createGeometryPyFromGeometryPtr(geometry *geometryPtr)
