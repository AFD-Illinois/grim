from gridHeaders cimport grid

cdef extern from "reconstruction.hpp" namespace "reconstruction":
  void reconstruct(const grid &prim,
              const int dir,
              grid &primLeft,
              grid &primRight,
              int &numReads,
              int &numWrites
             )
  void reconstructMM(const grid &prim,
                const int dir,
                grid &primLeft,
                grid &primRight,
                int &numReads,
                int &numWrites
               )
  void reconstructWENO5(const grid &prim,
                        const int dir,
                        grid &primLeft,
                        grid &primRight,
                        int &numReads,
                        int &numWrites
                       )
