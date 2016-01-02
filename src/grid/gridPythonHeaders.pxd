cdef extern from "grid.hpp":
  cdef cppclass grid:
    grid(int numVars)
    int numVars, numGhost, dim
    int iLocalStart, jLocalStart, kLocalStart
    int iLocalEnd,   jLocalEnd,   kLocalEnd
    int N1Local,     N2Local,     N3Local
    int N1Total,     N2Total,     N3Total
    int numGhostX1,  numGhostX2,  numGhostX3
    double dX1, dX2, dX3

    double *hostPtr
    void copyVarsToHostPtr()
    void communicate()
