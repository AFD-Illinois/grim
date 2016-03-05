"""
gridPy
======

Python interface to the grid class.
"""

import numpy as np
cimport numpy as np
from gridHeaders cimport grid
from gridHeaders cimport LOCATIONS_CENTER
from gridHeaders cimport LOCATIONS_LEFT
from gridHeaders cimport LOCATIONS_RIGHT
from gridHeaders cimport LOCATIONS_TOP
from gridHeaders cimport LOCATIONS_BOTTOM
from gridHeaders cimport LOCATIONS_FRONT
from gridHeaders cimport LOCATIONS_BACK

np.import_array()

# location macros
CENTER = LOCATIONS_CENTER
LEFT   = LOCATIONS_LEFT
RIGHT  = LOCATIONS_RIGHT
TOP    = LOCATIONS_TOP
BOTTOM = LOCATIONS_BOTTOM
FRONT  = LOCATIONS_FRONT
BACK   = LOCATIONS_BACK

cdef class gridPy(object):

  def __cinit__(self, int N1, int N2, int N3,
                      int dim,
                      int numVars,
                      int numGhost,
                      int periodicBoundariesX1,
                      int periodicBoundariesX2,
                      int periodicBoundariesX3
               ):
    self.gridPtr = new grid(N1, N2, N3, 
                            dim, numVars, numGhost, 
                            periodicBoundariesX1,
                            periodicBoundariesX2,
                            periodicBoundariesX3
                           )

  # Empty constructor needed to create gridPy from an existing gridPtr
  def __cinit__(self):
    self.gridPtr = NULL

  def __dealloc__(self):
    del self.gridPtr

  def __copy__(self):
    return gridPy(self.gridPtr.N1, 
                  self.gridPtr.N2,
                  self.gridPtr.N3,
                  self.gridPtr.dim, 
                  self.gridPtr.numVars,
                  self.gridPtr.numGhost,
                  self.gridPtr.periodicBoundariesX1,
                  self.gridPtr.periodicBoundariesX2,
                  self.gridPtr.periodicBoundariesX3
                 )

  def communicate(self):
    self.gridPtr.communicate()

  cdef grid* getGridPtr(self):
    return self.gridPtr

  cdef setGridPtr(self, grid *gridPtr):
    self.gridPtr = gridPtr  

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

  property N1Total:
    def __get__(self):
      return self.gridPtr.N1Total

  property N2Total:
    def __get__(self):
      return self.gridPtr.N2Total

  property N3Total:
    def __get__(self):
      return self.gridPtr.N3Total

  property numGhostX1:
    def __get__(self):
      return self.gridPtr.numGhostX1

  property numGhostX2:
    def __get__(self):
      return self.gridPtr.numGhostX2

  property numGhostX3:
    def __get__(self):
      return self.gridPtr.numGhostX3

  property dim:
    def __get__(self):
      return self.gridPtr.dim

  property shape:
    def __get__(self):
      return [self.gridPtr.numVars, \
              self.gridPtr.N3Total, \
              self.gridPtr.N2Total, \
              self.gridPtr.N1Total \
             ]

  def getVars(self):
    self.gridPtr.copyVarsToHostPtr()

    cdef int N1Total = self.gridPtr.N1Total
    cdef int N2Total = self.gridPtr.N2Total
    cdef int N3Total = self.gridPtr.N3Total
    cdef int numVars = self.gridPtr.numVars

    vars = \
      np.PyArray_SimpleNewFromData(4,
                                   [numVars, N3Total, N2Total, N1Total],
                                   np.NPY_DOUBLE,
                                   self.gridPtr.hostPtr
                                  )

    return vars

  def setVars(self, inputVars):
    cdef np.ndarray[double, ndim=4] vars = \
        np.ascontiguousarray(inputVars, dtype=np.float64)
    
    self.gridPtr.copyHostPtrToVars(&vars[0, 0, 0, 0])

cdef class coordinatesGridPy(object):

  def __cinit__(self, int N1, int N2, int N3,
                      int dim,
                      int numGhost,
                      double X1Start, double X1End,
                      double X2Start, double X2End,
                      double X3Start, double X3End
               ):
    self.coordGridPtr = new coordinatesGrid(N1, N2, N3, 
                                            dim, numGhost,
                                            X1Start, X1End,
                                            X2Start, X2End,
                                            X3Start, X3End
                                           )
  def __dealloc__(self):
    del self.coordGridPtr

  def __copy__(self):
    return coordinatesGridPy(self.coordGridPtr.N1,
                             self.coordGridPtr.N2,
                             self.coordGridPtr.N3,
                             self.coordGridPtr.dim,
                             self.coordGridPtr.numGhost,
                             self.coordGridPtr.X1Start, 
                             self.coordGridPtr.X1End,
                             self.coordGridPtr.X2Start,
                             self.coordGridPtr.X2End,
                             self.coordGridPtr.X3Start,
                             self.coordGridPtr.X3End 
                            )
  property N1:
    def __get__(self):
      return self.coordGridPtr.N1

  property N2:
    def __get__(self):
      return self.coordGridPtr.N2

  property N3:
    def __get__(self):
      return self.coordGridPtr.N3

  property numGhost:
    def __get__(self):
      return self.coordGridPtr.numGhost

  property dim:
    def __get__(self):
      return self.coordGridPtr.dim

  property dX1:
    def __get__(self):
      return self.coordGridPtr.dX1

  property dX2:
    def __get__(self):
      return self.coordGridPtr.dX2

  property dX3:
    def __get__(self):
      return self.coordGridPtr.dX3

  property X1Start:
    def __get__(self):
      return self.coordGridPtr.X1Start

  property X2Start:
    def __get__(self):
      return self.coordGridPtr.X2Start

  property X3Start:
    def __get__(self):
      return self.coordGridPtr.X3Start

  property X1End:
    def __get__(self):
      return self.coordGridPtr.X1End

  property X2End:
    def __get__(self):
      return self.coordGridPtr.X2End

  property X3End:
    def __get__(self):
      return self.coordGridPtr.X3End

  cdef coordinatesGrid* getGridPtr(self):
    return self.coordGridPtr

  def getCoords(self, int location):
    self.coordGridPtr.setXCoords(location)
    self.coordGridPtr.copyVarsToHostPtr()

    cdef int N1Total = self.coordGridPtr.N1Total
    cdef int N2Total = self.coordGridPtr.N2Total
    cdef int N3Total = self.coordGridPtr.N3Total

    cdef int numVars = 3

    vars = \
      np.PyArray_SimpleNewFromData(4,
                                   [numVars, N3Total, N2Total, N1Total],
                                   np.NPY_DOUBLE,
                                   self.coordGridPtr.hostPtr
                                  )

    return vars
