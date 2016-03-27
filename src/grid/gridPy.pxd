from gridHeaders cimport grid, coordinatesGrid

# Only the declarations here. Definitions in gridPy.pyx. This split is needed so
# that we can use share the gridPy module across different Cython modules using
# "from gridPy cimport gridPy".
# See "Sharing Extension Types" in
# http://docs.cython.org/src/userguide/sharing_declarations.html
cdef class gridPy(object):
  cdef grid *gridPtr
  cdef grid* getGridPtr(self)
  cdef setGridPtr(self, grid *gridPtr)
  cdef int usingExternalPtr

  @staticmethod
  cdef createGridPyFromGridPtr(grid *gridPtr)

cdef class coordinatesGridPy(object):
  cdef coordinatesGrid *coordGridPtr
  cdef coordinatesGrid* getGridPtr(self)
  cdef setGridPtr(self, coordinatesGrid *coordGridPtr)
