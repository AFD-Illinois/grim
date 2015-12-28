from gridPythonHeaders cimport grid

cdef class gridpy:
  cdef grid *gridPtr

  def __cinit__(self, int numVars):
    self.gridPtr = new grid(numVars)
    
  def __dealloc__(self):
    del self.gridPtr

  property numVars:
    def __get__(self):
      return self.gridPtr.numVars

#  property varsSoA:
#    def __get__(self):
#      return self.gridPtr.varsSoA
