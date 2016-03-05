from gridHeaders cimport grid, coordinatesGrid

cdef extern from "geometry.hpp":
  cdef enum:
    METRICS_MINKOWSKI             "metrics::MINKOWSKI"
    METRICS_MODIFIED_KERR_SCHILD  "metrics::MODIFIED_KERR_SCHILD"

  cdef cppclass geometry:
    geometry(const int metric,
             const double blackHoleSpin,
             const double hSlope,
             const coordinatesGrid &XCoordsGrid
            )
    int metric
    double blackHoleSpin
    double hSlope

    grid *gCovGrid
    grid *gConGrid
    grid *gGrid
    grid *alphaGrid
    grid *gammaUpDownDownGrid
    grid *xCoordsGrid

    void setgConGrid()
    void setgCovGrid()
    void setgGrid()
    void setalphaGrid()
    void setgammaUpDownDownGrid()
    void setxCoordsGrid()
