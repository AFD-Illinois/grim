import numpy as np
cimport numpy as np
from physicsHeaders cimport fluidElement
from gridPy cimport gridPy, coordinatesGridPy

cdef class fluidElementPy(object):
  cdef fluidElement *elemPtr
  cdef setElemPtr(self, fluidElement *elemPtr)
  cdef fluidElement* getElemPtr(self, fluidElement *elemPtr)
  cdef int usingExternalPtr

  @staticmethod
  cdef createFluidElementPyFromElemPtr(fluidElement *elemPtr)

