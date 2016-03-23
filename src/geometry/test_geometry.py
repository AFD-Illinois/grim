import mpi4py, petsc4py
from petsc4py import PETSc
import numpy as np
import pytest
import gridPy
import geometryPy

petsc4py.init()
petscComm  = petsc4py.PETSc.COMM_WORLD
comm = petscComm.tompi4py()
rank = comm.Get_rank()
numProcs = comm.Get_size()
PETSc.Sys.Print("Using %d procs" % numProcs)

N1  = int(pytest.config.getoption('N1'))
N2  = int(pytest.config.getoption('N2'))
N3  = int(pytest.config.getoption('N3'))
dim = int(pytest.config.getoption('dim'))

# Geometry parameters
blackHoleSpin = float(pytest.config.getoption('blackHoleSpin'))
hSlope        = float(pytest.config.getoption('hSlope'))
numGhost = 3
X1Start = 0.; X1End = 1.
X2Start = 0.; X2End = 1.
X3Start = 0.; X3End = 1.
periodicBoundariesX1 = False
periodicBoundariesX2 = False
periodicBoundariesX3 = False

XCoords = gridPy.coordinatesGridPy(N1, N2, N3,
                                   dim, numGhost,
                                   X1Start, X1End,
                                   X2Start, X2End,
                                   X3Start, X3End
                                  )
X1Coords, X2Coords, X3Coords = XCoords.getCoords(gridPy.CENTER)

geomMinkowski = geometryPy.geometryPy(geometryPy.MINKOWSKI,
                                      0., 0.,
                                      XCoords
                                     )

def test_minkowski_params():
  np.testing.assert_equal(N1,       geomMinkowski.N1)
  np.testing.assert_equal(N2,       geomMinkowski.N2)
  np.testing.assert_equal(N3,       geomMinkowski.N3)
  np.testing.assert_equal(dim,      geomMinkowski.dim)
  np.testing.assert_equal(numGhost, geomMinkowski.numGhost)

def test_minkowski_gCov():
  np.testing.assert_allclose(geomMinkowski.gCov[0][0], -1.)
  np.testing.assert_allclose(geomMinkowski.gCov[0][1],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[0][2],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[0][3],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[1][0],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[1][1],  1.)
  np.testing.assert_allclose(geomMinkowski.gCov[1][2],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[1][3],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[2][0],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[2][1],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[2][2],  1.)
  np.testing.assert_allclose(geomMinkowski.gCov[2][3],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[3][0],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[3][1],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[3][2],  0.)
  np.testing.assert_allclose(geomMinkowski.gCov[3][3],  1.)

def test_minkowski_gCon():
  np.testing.assert_allclose(geomMinkowski.gCon[0][0], -1.)
  np.testing.assert_allclose(geomMinkowski.gCon[0][1],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[0][2],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[0][3],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[1][0],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[1][1],  1.)
  np.testing.assert_allclose(geomMinkowski.gCon[1][2],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[1][3],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[2][0],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[2][1],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[2][2],  1.)
  np.testing.assert_allclose(geomMinkowski.gCon[2][3],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[3][0],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[3][1],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[3][2],  0.)
  np.testing.assert_allclose(geomMinkowski.gCon[3][3],  1.)

def test_minkowski_g():
  np.testing.assert_allclose(geomMinkowski.g, 1.)

def test_minkowski_alpha():
  np.testing.assert_allclose(geomMinkowski.g, 1.)

geomKerrSchild = geometryPy.geometryPy(geometryPy.MODIFIED_KERR_SCHILD,
                                       blackHoleSpin, hSlope,
                                       XCoords
                                      )
# From McKinney and Gammie, 2004
# Check if the coordinate transformations have been done correctly
r     = np.exp(X1Coords)
theta = np.pi*X2Coords + 0.5*(1. - hSlope)*np.sin(2.*np.pi*X2Coords)
phi   = 2*np.pi*X3Coords

sigma  = r**2. + (blackHoleSpin*np.cos(theta) )**2.
delta  = r**2. - 2*r + blackHoleSpin**2.
A      = (r**2. + blackHoleSpin**2.)**2.
sigmaMinus = r**2. - (blackHoleSpin*np.cos(theta) )**2.

# Coordinate transformation for log spacing in r and concentrating zones in the
# mid plane
dr_dX1       = np.exp(X1Coords)
dtheta_dX2   = np.pi*(1. + (1. - hSlope)*np.cos(2.*np.pi*X2Coords))
d2theta_dX22 = -2.*np.pi*np.pi*(1-hSlope)*np.sin(2.*np.pi*X2Coords);

N1Total = XCoords.N1Total
N2Total = XCoords.N2Total
N3Total = XCoords.N3Total
gCovCheck = np.zeros([4, 4, N3Total, N2Total, N1Total])
gConCheck = np.zeros([4, 4, N3Total, N2Total, N1Total])
gCheck    = np.zeros([N3Total, N2Total, N1Total])

gCovCheck[0][0] = -(1. - 2*r/sigma)    # dt^2
gCovCheck[0][1] = (2*r/sigma) * dr_dX1 # dt dX1
gCovCheck[0][2] = 0.                   # dt dX2
gCovCheck[0][3] = -(2.*blackHoleSpin*r*np.sin(theta)**2./sigma) # dt dphi
gCovCheck[1][0] = gCovCheck[0][1]
gCovCheck[1][1] = (1. + 2*r/sigma) * dr_dX1**2. # dX1 dX1
gCovCheck[1][2] = 0.
gCovCheck[1][3] = -blackHoleSpin * (1. + 2*r/sigma)*np.sin(theta)**2. \
                  * dr_dX1 # dX1 dphi
gCovCheck[2][0] = gCovCheck[0][2]
gCovCheck[2][1] = gCovCheck[1][2]
gCovCheck[2][2] = sigma * dtheta_dX2 * dtheta_dX2 # dX2 dX2
gCovCheck[2][3] = 0. # dX2 dphi
gCovCheck[3][0] = gCovCheck[0][3]
gCovCheck[3][1] = gCovCheck[1][3]
gCovCheck[3][2] = gCovCheck[2][3]
gCovCheck[3][3] =   np.sin(theta)**2. \
                  * (sigma +   blackHoleSpin**2. \
                             * (1. + 2.*r/sigma)*np.sin(theta)**2. \
                    ) # dphi dphi

gCovPerZone = np.zeros([4, 4])
for k in xrange(N3Total):
  for j in xrange(N2Total):
    for i in xrange(N1Total):
      gCovPerZone[0, 0] = gCovCheck[0][0][k, j, i]
      gCovPerZone[0, 1] = gCovCheck[0][1][k, j, i]
      gCovPerZone[0, 2] = gCovCheck[0][2][k, j, i]
      gCovPerZone[0, 3] = gCovCheck[0][3][k, j, i]
      gCovPerZone[1, 0] = gCovCheck[1][0][k, j, i]
      gCovPerZone[1, 1] = gCovCheck[1][1][k, j, i]
      gCovPerZone[1, 2] = gCovCheck[1][2][k, j, i]
      gCovPerZone[1, 3] = gCovCheck[1][3][k, j, i]
      gCovPerZone[2, 0] = gCovCheck[2][0][k, j, i]
      gCovPerZone[2, 1] = gCovCheck[2][1][k, j, i]
      gCovPerZone[2, 2] = gCovCheck[2][2][k, j, i]
      gCovPerZone[2, 3] = gCovCheck[2][3][k, j, i]
      gCovPerZone[3, 0] = gCovCheck[3][0][k, j, i]
      gCovPerZone[3, 1] = gCovCheck[3][1][k, j, i]
      gCovPerZone[3, 2] = gCovCheck[3][2][k, j, i]
      gCovPerZone[3, 3] = gCovCheck[3][3][k, j, i]

      gConPerZone     = np.linalg.inv(gCovPerZone)
      gCheck[k, j, i] = np.sqrt(-np.linalg.det(gCovPerZone))

      gConCheck[0][0][k, j, i] = gConPerZone[0, 0]
      gConCheck[0][1][k, j, i] = gConPerZone[0, 1]
      gConCheck[0][2][k, j, i] = gConPerZone[0, 2]
      gConCheck[0][3][k, j, i] = gConPerZone[0, 3]
      gConCheck[1][0][k, j, i] = gConPerZone[1, 0]
      gConCheck[1][1][k, j, i] = gConPerZone[1, 1]
      gConCheck[1][2][k, j, i] = gConPerZone[1, 2]
      gConCheck[1][3][k, j, i] = gConPerZone[1, 3]
      gConCheck[2][0][k, j, i] = gConPerZone[2, 0]
      gConCheck[2][1][k, j, i] = gConPerZone[2, 1]
      gConCheck[2][2][k, j, i] = gConPerZone[2, 2]
      gConCheck[2][3][k, j, i] = gConPerZone[2, 3]
      gConCheck[3][0][k, j, i] = gConPerZone[3, 0]
      gConCheck[3][1][k, j, i] = gConPerZone[3, 1]
      gConCheck[3][2][k, j, i] = gConPerZone[3, 2]
      gConCheck[3][3][k, j, i] = gConPerZone[3, 3]

alphaCheck = 1./np.sqrt(-gConCheck[0][0])

geomKerrSchild.computeConnectionCoeffs()

gammaUpDownDownCheck = np.zeros([4, 4, 4, N3Total, N2Total, N1Total])
gammaUpDownDownCheck[0][0][0] = 2.*r*sigmaMinus / sigma**3.
gammaUpDownDownCheck[0][0][1] = r * (2*r + sigma) * sigmaMinus / sigma**3.
gammaUpDownDownCheck[0][0][2] = -blackHoleSpin**2. * r * np.sin(2.*theta) \
                                * dtheta_dX2 / sigma**2.
gammaUpDownDownCheck[0][0][3] = -2. * blackHoleSpin * r * np.sin(theta)**2. \
                                * sigmaMinus / sigma**3.

gammaUpDownDownCheck[0][1][0] = gammaUpDownDownCheck[0][0][1]
gammaUpDownDownCheck[0][1][1] = 2.*r**2.*(r**4. + r*sigmaMinus
                                          - (blackHoleSpin*np.cos(theta))**4.
                                         ) / sigma**3.
gammaUpDownDownCheck[0][1][2] = -blackHoleSpin**2. * r**2. * np.sin(2.*theta) \
                                * dtheta_dX2 / sigma**2.
gammaUpDownDownCheck[0][1][3] = blackHoleSpin * r * (-r*(r**3. + 2*sigmaMinus)
                                                     + ( blackHoleSpin 
                                                         * np.cos(theta)
                                                       )**4.
                                                    ) * np.sin(theta)**2. \
                                                    / sigma**3.
                              
gammaUpDownDownCheck[0][2][0] = gammaUpDownDownCheck[0][0][2]
gammaUpDownDownCheck[0][2][1] = gammaUpDownDownCheck[0][1][2]
gammaUpDownDownCheck[0][2][2] = -2. * r**2. * dtheta_dX2**2. / sigma
gammaUpDownDownCheck[0][2][3] = blackHoleSpin**3. * r * np.sin(theta)**2. \
                                * np.sin(2.*theta) * dtheta_dX2 / sigma**2.

gammaUpDownDownCheck[0][3][0] = gammaUpDownDownCheck[0][0][3]
gammaUpDownDownCheck[0][3][1] = gammaUpDownDownCheck[0][1][3]
gammaUpDownDownCheck[0][3][2] = gammaUpDownDownCheck[0][2][3]
gammaUpDownDownCheck[0][3][3] = 2.*r*np.sin(theta)**2. \
                                * (-r*sigma**2. +
                                   blackHoleSpin**2.*np.sin(theta)**2.*sigmaMinus
                                  ) / sigma**3.

gammaUpDownDownCheck[1][0][0] = (blackHoleSpin**2. + r*(-2. + r)) \
                                * sigmaMinus / (r * sigma**3.)
gammaUpDownDownCheck[1][0][1] = sigmaMinus \
                                * ( -2.*r + (blackHoleSpin*np.sin(theta))**2.) \
                                / sigma**3.
gammaUpDownDownCheck[1][0][2] = 0.
gammaUpDownDownCheck[1][0][3] = -blackHoleSpin * np.sin(theta)**2. \
                                * (blackHoleSpin**2. + r*(-2. + r)) * sigmaMinus \
                                / (r * sigma**3.)

gammaUpDownDownCheck[1][1][0] = gammaUpDownDownCheck[1][0][1]
gammaUpDownDownCheck[1][1][1] = \
  (  r**4.*(-2. + r)*(1. + r) 
   + blackHoleSpin**2. * (  blackHoleSpin**2.*r*(1. + 3.*r)*np.cos(theta)**4. \
                          + (blackHoleSpin*np.cos(theta))**4. * np.cos(theta)**2. \
                          + r**3.*np.sin(theta)**2. \
                          + r*np.cos(theta)**2. \
                            *(2.*r + 3.*r**3. - (blackHoleSpin*np.sin(theta))**2.)
                         )
  ) / sigma**3.
gammaUpDownDownCheck[1][1][2] = -blackHoleSpin**2. * dtheta_dX2 \
                               * np.sin(2.*theta) \
                               / (blackHoleSpin**2. + 2.*r**2.
                                  + blackHoleSpin**2.*np.cos(2.*theta) 
                                 )
gammaUpDownDownCheck[1][1][3] = \
  blackHoleSpin * np.sin(theta)**2. * (blackHoleSpin**4. * r * np.cos(theta)**4.
                                       + r**2*(2.*r + r**3. 
                                               -(blackHoleSpin*np.sin(theta))**2.
                                              ) \
                                       + (blackHoleSpin*np.cos(theta))**2. \
                                       * (2.*r*(-1. + r**2.)
                                          + (blackHoleSpin*np.sin(theta))**2.
                                         )
                                      ) / sigma**3.

gammaUpDownDownCheck[1][2][0] = gammaUpDownDownCheck[1][0][2]
gammaUpDownDownCheck[1][2][1] = gammaUpDownDownCheck[1][1][2]
gammaUpDownDownCheck[1][2][2] = -(blackHoleSpin**2. + r*(-2. + r)) \
                               * dtheta_dX2**2. / sigma
gammaUpDownDownCheck[1][2][3] = 0.

gammaUpDownDownCheck[1][3][0] = gammaUpDownDownCheck[1][0][3]
gammaUpDownDownCheck[1][3][1] = gammaUpDownDownCheck[1][1][3]
gammaUpDownDownCheck[1][3][2] = gammaUpDownDownCheck[1][2][3]
gammaUpDownDownCheck[1][3][3] = \
  -(blackHoleSpin**2. + r*(-2. + r) ) * np.sin(theta)**2. \
                                * (r * sigma**2. -
                                   blackHoleSpin**2.*sigmaMinus*np.sin(theta)**2.
                                  ) / (r * sigma**3.)

gammaUpDownDownCheck[2][0][0] = -blackHoleSpin**2. * r * np.sin(2.*theta) \
                                / sigma**3. / dtheta_dX2
gammaUpDownDownCheck[2][0][1] = r * gammaUpDownDownCheck[2][0][0]
gammaUpDownDownCheck[2][0][2] = 0.
gammaUpDownDownCheck[2][0][3] = blackHoleSpin*r*(blackHoleSpin**2. + r**2.) \
                                * np.sin(2.*theta) / sigma**3. / dtheta_dX2

gammaUpDownDownCheck[2][1][0] = gammaUpDownDownCheck[2][0][1]
gammaUpDownDownCheck[2][1][1] = r**2. * gammaUpDownDownCheck[2][0][0]
gammaUpDownDownCheck[2][1][2] = r**2. / sigma
gammaUpDownDownCheck[2][1][3] = (blackHoleSpin*r*np.cos(theta)*np.sin(theta)
                                 *(r**3.*(2. + r) 
                                    + blackHoleSpin**2.
                                    *( 2.*r*(1. + r)*np.cos(theta)**2.
                                      + blackHoleSpin**2.*np.cos(theta)**4. 
                                      + 2.*r*np.sin(theta)**2.
                                     )
                                  )
                                ) / sigma**3. / dtheta_dX2

gammaUpDownDownCheck[2][2][0] = gammaUpDownDownCheck[2][0][2]
gammaUpDownDownCheck[2][2][1] = gammaUpDownDownCheck[2][1][2]
gammaUpDownDownCheck[2][2][2] = -blackHoleSpin**2.*np.cos(theta)*np.sin(theta) \
                                *dtheta_dX2/sigma + d2theta_dX22/dtheta_dX2
gammaUpDownDownCheck[2][2][3] = 0.

gammaUpDownDownCheck[2][3][0] = gammaUpDownDownCheck[2][0][3]
gammaUpDownDownCheck[2][3][1] = gammaUpDownDownCheck[2][1][3]
gammaUpDownDownCheck[2][3][2] = gammaUpDownDownCheck[2][2][3]
gammaUpDownDownCheck[2][3][3] = \
  -np.cos(theta)*np.sin(theta) \
  *(sigma**3. + (blackHoleSpin*np.sin(theta))**2. \
                * sigma*(r*(4. + r) + (blackHoleSpin*np.cos(theta)**2.)) \
              + 2.*r*(blackHoleSpin * np.sin(theta))**4. \
   ) / sigma**3. / dtheta_dX2

gammaUpDownDownCheck[3][0][0] = blackHoleSpin * sigmaMinus / sigma**3.
gammaUpDownDownCheck[3][0][1] = r * gammaUpDownDownCheck[3][0][0]
gammaUpDownDownCheck[3][0][2] = -2.*blackHoleSpin*r*np.cos(theta) \
                                * dtheta_dX2 / (np.sin(theta) * sigma**2.)
gammaUpDownDownCheck[3][0][3] = -blackHoleSpin**2. * np.sin(theta)**2. \
                                * sigmaMinus / sigma**3.

gammaUpDownDownCheck[3][1][0] = gammaUpDownDownCheck[3][0][1]
gammaUpDownDownCheck[3][1][1] = blackHoleSpin * r**2. * sigmaMinus \
                                / sigma**3.
gammaUpDownDownCheck[3][1][2] = -2.*blackHoleSpin*r \
                                *(blackHoleSpin**2. + 2.*r*(2. + r)
                                  + blackHoleSpin**2. * np.cos(2.*theta)
                                 ) * np.cos(theta) * dtheta_dX2 \
                                / (np.sin(theta) \
                                   * (blackHoleSpin**2. + 2.*r**2.
                                      + blackHoleSpin**2.*np.cos(2.*theta)
                                     )**2.
                                  )
gammaUpDownDownCheck[3][1][3] = \
  r*(r*sigma**2. - (blackHoleSpin*np.sin(theta))**2.*sigmaMinus)/sigma**3.

gammaUpDownDownCheck[3][2][0] = gammaUpDownDownCheck[3][0][2]
gammaUpDownDownCheck[3][2][1] = gammaUpDownDownCheck[3][1][2]
gammaUpDownDownCheck[3][2][2] = -blackHoleSpin * r * dtheta_dX2**2./sigma
gammaUpDownDownCheck[3][2][3] = \
    dtheta_dX2*(.25*(blackHoleSpin**2.
                     + 2.*r**2. + blackHoleSpin**2.*np.cos(2.*theta)
                    )**2. * np.cos(theta)/np.sin(theta)
                + blackHoleSpin**2. * r * np.sin(2.*theta)
               )/sigma**2.

gammaUpDownDownCheck[3][3][0] = gammaUpDownDownCheck[3][0][3]
gammaUpDownDownCheck[3][3][1] = gammaUpDownDownCheck[3][1][3]
gammaUpDownDownCheck[3][3][2] = gammaUpDownDownCheck[3][2][3]
gammaUpDownDownCheck[3][3][3] = \
    (-blackHoleSpin * r * np.sin(theta)**2. * sigma**2. \
     + blackHoleSpin**3. * np.sin(theta)**4. * sigmaMinus) / sigma**3.

def test_modifiedKerrSchild_params():                       
  np.testing.assert_equal(N1,       geomKerrSchild.N1)      
  np.testing.assert_equal(N2,       geomKerrSchild.N2)
  np.testing.assert_equal(N3,       geomKerrSchild.N3)
  np.testing.assert_equal(dim,      geomKerrSchild.dim)
  np.testing.assert_equal(numGhost, geomKerrSchild.numGhost)

def test_modifiedKerrSchild_xCoords():

  np.testing.assert_allclose(r,     geomKerrSchild.xCoords[0])
  np.testing.assert_allclose(theta, geomKerrSchild.xCoords[1])
  np.testing.assert_allclose(phi,   geomKerrSchild.xCoords[2])

def test_modifiedKerrSchild_gCov():
  np.testing.assert_allclose(gCovCheck[0][0], geomKerrSchild.gCov[0][0])
  np.testing.assert_allclose(gCovCheck[0][1], geomKerrSchild.gCov[0][1])
  np.testing.assert_allclose(gCovCheck[0][2], geomKerrSchild.gCov[0][2])
  np.testing.assert_allclose(gCovCheck[0][3], geomKerrSchild.gCov[0][3])
  np.testing.assert_allclose(gCovCheck[1][0], geomKerrSchild.gCov[1][0])
  np.testing.assert_allclose(gCovCheck[1][1], geomKerrSchild.gCov[1][1])
  np.testing.assert_allclose(gCovCheck[1][2], geomKerrSchild.gCov[1][2])
  np.testing.assert_allclose(gCovCheck[1][3], geomKerrSchild.gCov[1][3])
  np.testing.assert_allclose(gCovCheck[2][0], geomKerrSchild.gCov[2][0])
  np.testing.assert_allclose(gCovCheck[2][1], geomKerrSchild.gCov[2][1])
  np.testing.assert_allclose(gCovCheck[2][2], geomKerrSchild.gCov[2][2])
  np.testing.assert_allclose(gCovCheck[2][3], geomKerrSchild.gCov[2][3])
  np.testing.assert_allclose(gCovCheck[3][0], geomKerrSchild.gCov[3][0])
  np.testing.assert_allclose(gCovCheck[3][1], geomKerrSchild.gCov[3][1])
  np.testing.assert_allclose(gCovCheck[3][2], geomKerrSchild.gCov[3][2])
  np.testing.assert_allclose(gCovCheck[3][3], geomKerrSchild.gCov[3][3])

def test_modifiedKerrSchild_gCon():
  np.testing.assert_allclose(gConCheck[0][0], geomKerrSchild.gCon[0][0])
  np.testing.assert_allclose(gConCheck[0][1], geomKerrSchild.gCon[0][1])
  np.testing.assert_allclose(gConCheck[0][2], geomKerrSchild.gCon[0][2])
  np.testing.assert_allclose(gConCheck[0][3], geomKerrSchild.gCon[0][3],
                             atol=1e-14
                            )
  np.testing.assert_allclose(gConCheck[1][0], geomKerrSchild.gCon[1][0])
  np.testing.assert_allclose(gConCheck[1][1], geomKerrSchild.gCon[1][1])
  np.testing.assert_allclose(gConCheck[1][2], geomKerrSchild.gCon[1][2])
  np.testing.assert_allclose(gConCheck[1][3], geomKerrSchild.gCon[1][3])
  np.testing.assert_allclose(gConCheck[2][0], geomKerrSchild.gCon[2][0])
  np.testing.assert_allclose(gConCheck[2][1], geomKerrSchild.gCon[2][1])
  np.testing.assert_allclose(gConCheck[2][2], geomKerrSchild.gCon[2][2])
  np.testing.assert_allclose(gConCheck[2][3], geomKerrSchild.gCon[2][3])
  np.testing.assert_allclose(gConCheck[3][0], geomKerrSchild.gCon[3][0], 
                             atol=1e-14
                            )
  np.testing.assert_allclose(gConCheck[3][1], geomKerrSchild.gCon[3][1])
  np.testing.assert_allclose(gConCheck[3][2], geomKerrSchild.gCon[3][2])
  np.testing.assert_allclose(gConCheck[3][3], geomKerrSchild.gCon[3][3])
  
def test_modifiedKerrSchild_g():
  np.testing.assert_allclose(gCheck, geomKerrSchild.g)

def test_modifiedKerrSchild_alpha():
  np.testing.assert_allclose(alphaCheck, geomKerrSchild.alpha)

def test_modifiedKerrSchild_gammaUpDownDown():
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][0][0],
                             geomKerrSchild.gammaUpDownDown[0][0][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][0][1],
                             geomKerrSchild.gammaUpDownDown[0][0][1]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][0][2],
                             geomKerrSchild.gammaUpDownDown[0][0][2],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][0][3],
                             geomKerrSchild.gammaUpDownDown[0][0][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][1][0],
                             geomKerrSchild.gammaUpDownDown[0][1][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][1][1],
                             geomKerrSchild.gammaUpDownDown[0][1][1]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][1][2],
                             geomKerrSchild.gammaUpDownDown[0][1][2],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][1][3],
                             geomKerrSchild.gammaUpDownDown[0][1][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][2][0],
                             geomKerrSchild.gammaUpDownDown[0][2][0],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][2][1],
                             geomKerrSchild.gammaUpDownDown[0][2][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][2][2],
                             geomKerrSchild.gammaUpDownDown[0][2][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][2][3],
                             geomKerrSchild.gammaUpDownDown[0][2][3],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][3][0],
                             geomKerrSchild.gammaUpDownDown[0][3][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][3][1],
                             geomKerrSchild.gammaUpDownDown[0][3][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][3][2],
                             geomKerrSchild.gammaUpDownDown[0][3][2],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[0][3][3],
                             geomKerrSchild.gammaUpDownDown[0][3][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][0][0],
                             geomKerrSchild.gammaUpDownDown[1][0][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][0][1],
                             geomKerrSchild.gammaUpDownDown[1][0][1]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][0][2],
                             geomKerrSchild.gammaUpDownDown[1][0][2],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][0][3],
                             geomKerrSchild.gammaUpDownDown[1][0][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][1][0],
                             geomKerrSchild.gammaUpDownDown[1][1][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][1][1],
                             geomKerrSchild.gammaUpDownDown[1][1][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][1][2],
                             geomKerrSchild.gammaUpDownDown[1][1][2],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][1][3],
                             geomKerrSchild.gammaUpDownDown[1][1][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][2][0],
                             geomKerrSchild.gammaUpDownDown[1][2][0],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][2][1],
                             geomKerrSchild.gammaUpDownDown[1][2][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][2][2],
                             geomKerrSchild.gammaUpDownDown[1][2][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][2][3],
                             geomKerrSchild.gammaUpDownDown[1][2][3],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][3][0],
                             geomKerrSchild.gammaUpDownDown[1][3][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][3][1],
                             geomKerrSchild.gammaUpDownDown[1][3][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][3][2],
                             geomKerrSchild.gammaUpDownDown[1][3][2],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[1][3][3],
                             geomKerrSchild.gammaUpDownDown[1][3][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][0][0],
                             geomKerrSchild.gammaUpDownDown[2][0][0],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][0][1],
                             geomKerrSchild.gammaUpDownDown[2][0][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][0][2],
                             geomKerrSchild.gammaUpDownDown[2][0][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][0][3],
                             geomKerrSchild.gammaUpDownDown[2][0][3],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][1][0],
                             geomKerrSchild.gammaUpDownDown[2][1][0],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][1][1],
                             geomKerrSchild.gammaUpDownDown[2][1][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][1][2],
                             geomKerrSchild.gammaUpDownDown[2][1][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][1][3],
                             geomKerrSchild.gammaUpDownDown[2][1][3],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][2][0],
                             geomKerrSchild.gammaUpDownDown[2][2][0],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][2][1],
                             geomKerrSchild.gammaUpDownDown[2][2][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][2][2],
                             geomKerrSchild.gammaUpDownDown[2][2][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][2][3],
                             geomKerrSchild.gammaUpDownDown[2][2][3],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][3][0],
                             geomKerrSchild.gammaUpDownDown[2][3][0],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][3][1],
                             geomKerrSchild.gammaUpDownDown[2][3][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][3][2],
                             geomKerrSchild.gammaUpDownDown[2][3][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[2][3][3],
                             geomKerrSchild.gammaUpDownDown[2][3][3],
                             atol=3e-3
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][0][0],
                             geomKerrSchild.gammaUpDownDown[3][0][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][0][1],
                             geomKerrSchild.gammaUpDownDown[3][0][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][0][2],
                             geomKerrSchild.gammaUpDownDown[3][0][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][0][3],
                             geomKerrSchild.gammaUpDownDown[3][0][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][1][0],
                             geomKerrSchild.gammaUpDownDown[3][1][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][1][1],
                             geomKerrSchild.gammaUpDownDown[3][1][1],
                             atol=1e-7
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][1][2],
                             geomKerrSchild.gammaUpDownDown[3][1][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][1][3],
                             geomKerrSchild.gammaUpDownDown[3][1][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][2][0],
                             geomKerrSchild.gammaUpDownDown[3][2][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][2][1],
                             geomKerrSchild.gammaUpDownDown[3][2][1]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][2][2],
                             geomKerrSchild.gammaUpDownDown[3][2][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][2][3],
                             geomKerrSchild.gammaUpDownDown[3][2][3]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][3][0],
                             geomKerrSchild.gammaUpDownDown[3][3][0]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][3][1],
                             geomKerrSchild.gammaUpDownDown[3][3][1]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][3][2],
                             geomKerrSchild.gammaUpDownDown[3][3][2]
                            )
  np.testing.assert_allclose(          gammaUpDownDownCheck[3][3][3],
                             geomKerrSchild.gammaUpDownDown[3][3][3]
                            )
