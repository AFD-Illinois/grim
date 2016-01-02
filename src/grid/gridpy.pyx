import numpy as np
cimport numpy as np
from gridPythonHeaders cimport grid

np.import_array()

cdef class gridpy(object):
  cdef grid *gridPtr
  cdef vars

  def __cinit__(self, int numVars):
    self.gridPtr = new grid(numVars)
    
  def __dealloc__(self):
    del self.gridPtr

  def __copy__(self):
    return gridpy(self.gridPtr.numVars)

  def communicate(self):
    self.gridPtr.communicate()

  property numVars:
    def __get__(self):
      return self.gridPtr.numVars

  property iLocalStart:
    def __get__(self):
      return self.gridPtr.iLocalStart

  property jLocalStart:
    def __get__(self):
      return self.gridPtr.jLocalStart

  property kLocalStart:
    def __get__(self):
      return self.gridPtr.kLocalStart

  property iLocalEnd:
    def __get__(self):
      return self.gridPtr.iLocalEnd

  property jLocalEnd:
    def __get__(self):
      return self.gridPtr.jLocalEnd

  property kLocalEnd:
    def __get__(self):
      return self.gridPtr.kLocalEnd

  property N1Local:
    def __get__(self):
      return self.gridPtr.N1Local

  property N2Local:
    def __get__(self):
      return self.gridPtr.N2Local

  property N3Local:
    def __get__(self):
      return self.gridPtr.N3Local

  property numGhostX1:
    def __get__(self):
      return self.gridPtr.numGhostX1

  property numGhostX2:
    def __get__(self):
      return self.gridPtr.numGhostX2

  property numGhostX3:
    def __get__(self):
      return self.gridPtr.numGhostX3

  property dX1:
    def __get__(self):
      return self.gridPtr.dX1

  property dX2:
    def __get__(self):
      return self.gridPtr.dX2

  property dX3:
    def __get__(self):
      return self.gridPtr.dX3

  property dim:
    def __get__(self):
      return self.gridPtr.dim

  def getVars(self):
      self.gridPtr.copyVarsToHostPtr()

      cdef int N1Total = self.gridPtr.N1Total
      cdef int N2Total = self.gridPtr.N2Total
      cdef int N3Total = self.gridPtr.N3Total
      cdef int numVars = self.gridPtr.numVars

      self.vars = \
        np.PyArray_SimpleNewFromData(4,
                                     [numVars, N3Total, N2Total, N1Total],
                                     np.NPY_DOUBLE,
                                     self.gridPtr.hostPtr
                                    )

      return self.vars

  def setVars(self, inputVars):
    cdef np.ndarray[double, ndim=4] vars = \
        np.ascontiguousarray(inputVars, dtype=np.float64)
    
    self.gridPtr.copyHostPtrToVars(&vars[0, 0, 0, 0])

