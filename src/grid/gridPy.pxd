from gridPythonHeaders cimport grid

cdef class gridPy(object):
  cdef grid *gridPtr
  cdef grid* getGridPtr(self)
