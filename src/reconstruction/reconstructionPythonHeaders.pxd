from gridPythonHeaders cimport grid

cdef extern from "reconstruction.hpp" namespace "reconstruction":
  void reconstruct(const grid &prim,
                   const int dir,
                   grid &primLeft,
                   grid &primRight
                  )
  void reconstructMM(const grid &prim,
                     const int dir,
                     grid &primLeft,
                     grid &primRight
                    )
  void reconstructWENO5(const grid &prim,
                        const int dir,
                        grid &primLeft,
                        grid &primRight
                       )
