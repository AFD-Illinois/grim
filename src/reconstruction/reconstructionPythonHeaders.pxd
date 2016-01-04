from gridPythonHeaders cimport grid

cdef extern from "reconstruction.hpp" namespace "reconstruction":
  void reconstruct(const grid &prim,
                   const int dir,
                   grid &primLeft,
                   grid &primRight
                  )
