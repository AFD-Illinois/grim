cdef extern from "grid.hpp":
  cdef cppclass grid:
    grid(int numVars)
    int numVars
