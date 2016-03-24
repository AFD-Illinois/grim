import numpy as np
cimport numpy as np
from gridPy cimport gridPy, coordinatesGridPy
from geometryHeaders cimport geometry
from geometryHeaders cimport METRICS_MINKOWSKI
from geometryHeaders cimport METRICS_MODIFIED_KERR_SCHILD

np.import_array()

# Metric macros
MINKOWSKI             = METRICS_MINKOWSKI
MODIFIED_KERR_SCHILD  = METRICS_MODIFIED_KERR_SCHILD

cdef class geometryPy(object):

  def __cinit__(self, int metric, 
                      double blackHoleSpin, 
                      double hSlope,
                      coordinatesGridPy XCoordsGrid
               ):
    self.geometryPtr = new geometry(metric,
                                    blackHoleSpin,
                                    hSlope,
                                    XCoordsGrid.getGridPtr()[0]
                                   )
    
    self.geometryPtr.setgCovGrid()
    self.geometryPtr.setgConGrid()
    self.geometryPtr.setgGrid()
    self.geometryPtr.setalphaGrid()
    self.geometryPtr.setxCoordsGrid()

    # TODO: make a function createGridPyFromGridPtr that will directly return a
    # gridPy class given a grid *gridPtr
    gCovGridPy = gridPy()
    gCovGridPy.setGridPtr(self.geometryPtr.gCovGrid)

    gConGridPy = gridPy()
    gConGridPy.setGridPtr(self.geometryPtr.gConGrid)

    gGridPy = gridPy()
    gGridPy.setGridPtr(self.geometryPtr.gGrid)

    alphaGridPy = gridPy()
    alphaGridPy.setGridPtr(self.geometryPtr.alphaGrid)

    xCoordsGridPy = gridPy()
    xCoordsGridPy.setGridPtr(self.geometryPtr.xCoordsGrid)

    self.gCov = gCovGridPy.getVars()
    self.gCov = self.gCov.reshape([4, 4, 
                                   self.gCov.shape[1],
                                   self.gCov.shape[2],
                                   self.gCov.shape[3]
                                  ]
                                 )
                             
    self.gCon = gConGridPy.getVars()
    self.gCon = self.gCon.reshape([4, 4, 
                                   self.gCon.shape[1],
                                   self.gCon.shape[2],
                                   self.gCon.shape[3]
                                  ]
                                 )

    self.g = gGridPy.getVars()
    self.g = self.g.reshape([self.g.shape[1],
                             self.g.shape[2],
                             self.g.shape[3]
                            ]
                           )
    self.alpha = alphaGridPy.getVars()
    self.alpha = self.alpha.reshape([self.alpha.shape[1],
                                     self.alpha.shape[2],
                                     self.alpha.shape[3]
                                    ]
                                   )
    self.xCoords = xCoordsGridPy.getVars()

  def __dealloc__(self):
    del self.geometryPtr

  def computeConnectionCoeffs(self):
    self.geometryPtr.computeConnectionCoeffs()
    self.geometryPtr.setgammaUpDownDownGrid()
    gammaUpDownDownGridPy = gridPy()
    gammaUpDownDownGridPy.setGridPtr(self.geometryPtr.gammaUpDownDownGrid)
    self.gammaUpDownDown = gammaUpDownDownGridPy.getVars()
    self.gammaUpDownDown = \
          self.gammaUpDownDown.reshape([4, 4, 4, 
                                        self.gammaUpDownDown.shape[1],
                                        self.gammaUpDownDown.shape[2],
                                        self.gammaUpDownDown.shape[3]
                                       ]
                                      )

  property N1:
    def __get__(self):
      return self.geometryPtr.N1

  property N2:
    def __get__(self):
      return self.geometryPtr.N2

  property N3:
    def __get__(self):
      return self.geometryPtr.N3

  property dim:
    def __get__(self):
      return self.geometryPtr.dim

  property numGhost:
    def __get__(self):
      return self.geometryPtr.numGhost

  property gCon:
    def __get__(self):
      return self.gCon

  property gCov:
    def __get__(self):
      return self.gCov

  property g:
    def __get__(self):
      return self.g

  property alpha:
    def __get__(self):
      return self.alpha

  property gammaUpDownDown:
    def __get__(self):
      return self.gammaUpDownDown
  
  property xCoords:
    def __get__(self):
      return self.xCoords

  property metric:
    def __get__(self):
      return self.geometryPtr.metric

  property blackHoleSpin:
    def __get__(self):
      return self.geometryPtr.blackHoleSpin

  property hSlope:
    def __get__(self):
      return self.geometryPtr.hSlope
