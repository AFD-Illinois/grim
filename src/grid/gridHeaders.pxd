cdef extern from "grid.hpp":
  cdef enum:
      LOCATIONS_BACK   "locations::BACK"
      LOCATIONS_CENTER "locations::CENTER"
      LOCATIONS_LEFT   "locations::LEFT"
      LOCATIONS_RIGHT  "locations::RIGHT"
      LOCATIONS_TOP    "locations::TOP"
      LOCATIONS_BOTTOM "locations::BOTTOM"
      LOCATIONS_FRONT  "locations::FRONT"

cdef extern from "grid.hpp":
  cdef cppclass grid:
    grid(int N1, int N2, int N3, 
         int dim, 
         int numVars,
         int numGhost,
         int periodicBoundariesX1,
         int periodicBoundariesX2,
         int periodicBoundariesX3
        )
    int N1, N2, N3
    int numVars, numGhost, dim
    int iLocalStart, jLocalStart, kLocalStart
    int iLocalEnd,   jLocalEnd,   kLocalEnd
    int N1Local,     N2Local,     N3Local
    int N1Total,     N2Total,     N3Total
    int numGhostX1,  numGhostX2,  numGhostX3
    double boundaryLeft,    boundaryRight
    double boundaryTop,     boundaryBottom
    double boundaryFront,   boundaryBack
    int periodicBoundariesX1
    int periodicBoundariesX2
    int periodicBoundariesX3

    double *hostPtr
    void communicate()
    void copyVarsToHostPtr()
    void copyHostPtrToVars(const double *hostPtr)

  cdef cppclass coordinatesGrid:
    coordinatesGrid(int N1, int N2, int N3, 
                    int dim, 
                    int numGhost,
                    double X1Start, double X1End,
                    double X2Start, double X2End,
                    double X3Start, double X3End
                   )
    int N1, N2, N3
    int N1Total, N2Total, N3Total
    int numGhost, dim
    double dX1, dX2, dX3
    double X1Start, X1End
    double X2Start, X2End
    double X3Start, X3End

    double *hostPtr
    void setXCoords(const int location)
    void copyVarsToHostPtr()
