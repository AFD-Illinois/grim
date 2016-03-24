import mpi4py, petsc4py
import numpy as np
import pytest
import gridPy

petsc4py.init()
petscComm  = petsc4py.PETSc.COMM_WORLD
comm = petscComm.tompi4py()
rank = comm.Get_rank()
numProcs = comm.Get_size()

N1  = int(pytest.config.getoption('N1'))
N2  = int(pytest.config.getoption('N2'))
N3  = int(pytest.config.getoption('N3'))
dim = int(pytest.config.getoption('dim'))
numVars = 1
numGhost = 3
X1Start = 0.; X1End = 1.
X2Start = 0.; X2End = 1.
X3Start = 0.; X3End = 1.
periodicBoundariesX1 = True
periodicBoundariesX2 = True
periodicBoundariesX3 = True

prim = gridPy.gridPy(N1, N2, N3, 
                     dim, numVars, numGhost,
                     periodicBoundariesX1,
                     periodicBoundariesX2,
                     periodicBoundariesX3
                    )

XCoords = gridPy.coordinatesGridPy(N1, N2, N3,
                                   dim, numGhost,
                                   X1Start, X1End,
                                   X2Start, X2End,
                                   X3Start, X3End
                                  )
X1Coords, X2Coords, X3Coords = XCoords.getCoords(gridPy.CENTER)

numGhostX1 = prim.numGhostX1
numGhostX2 = prim.numGhostX2
numGhostX3 = prim.numGhostX3

if (dim==1):
  domainX1Start =  numGhostX1
  domainX1End   = -numGhostX1

  domainX2Start =  0
  domainX2End   =  1

  domainX3Start =  0
  domainX3End   =  1

elif (dim==2):

  domainX1Start =  numGhostX1
  domainX1End   = -numGhostX1

  domainX2Start =  numGhostX2
  domainX2End   = -numGhostX2

  domainX3Start =  0
  domainX3End   =  1

elif (dim==3):

  domainX1Start =  numGhostX1
  domainX1End   = -numGhostX1

  domainX2Start =  numGhostX2
  domainX2End   = -numGhostX2

  domainX3Start =  numGhostX3
  domainX3End   = -numGhostX3

X1CoordsBulk = X1Coords[domainX3Start:domainX3End, \
                        domainX2Start:domainX2End, \
                        domainX1Start:domainX1End \
                       ]
X2CoordsBulk = X2Coords[domainX3Start:domainX3End, \
                        domainX2Start:domainX2End, \
                        domainX1Start:domainX1End \
                       ]
X3CoordsBulk = X3Coords[domainX3Start:domainX3End, \
                        domainX2Start:domainX2End, \
                        domainX1Start:domainX1End \
                       ]

varsNumpy = prim.getVars()
for var in xrange(prim.numVars):
  varsNumpy[var, \
            domainX3Start:domainX3End, \
            domainX2Start:domainX2End, \
            domainX1Start:domainX1End  \
           ] = \
    np.sin(2.*np.pi*X1CoordsBulk) \
  * np.sin(2.*np.pi*X2CoordsBulk) \
  * np.sin(2.*np.pi*X3CoordsBulk)

prim.setVars(varsNumpy)
prim.communicate()

varsNumpy = prim.getVars()
if (dim==1):
  d_dX1 = np.gradient(varsNumpy[0, 0, 0, :], XCoords.dX1)
  d_dX2 = 0.
  d_dX3 = 0.

elif (dim==2):
  d_dX2, d_dX1 = np.gradient(varsNumpy[0, 0, :, :], \
                             XCoords.dX2, XCoords.dX1 \
                            )
  d_dX3 = 0.
elif (dim==3):
  d_dX3, d_dX2, d_dX1  = \
      np.gradient(varsNumpy[0, :, :, :], \
                  XCoords.dX3, XCoords.dX2, XCoords.dX1 \
                 )

d_dX1Analytical =   2.*np.pi * np.cos(2.*np.pi*X1Coords) \
                             * np.sin(2.*np.pi*X2Coords) \
                             * np.sin(2.*np.pi*X3Coords)

d_dX2Analytical =              np.sin(2.*np.pi*X1Coords) \
                  * 2.*np.pi * np.cos(2.*np.pi*X2Coords) \
                             * np.sin(2.*np.pi*X3Coords)

d_dX3Analytical =              np.sin(2.*np.pi*X1Coords) \
                             * np.sin(2.*np.pi*X2Coords) \
                  * 2.*np.pi * np.cos(2.*np.pi*X3Coords)

errorX1 = np.abs((d_dX1 - d_dX1Analytical)[domainX3Start:domainX3End, \
                                           domainX2Start:domainX2End, \
                                           domainX1Start:domainX1End])

errorX2 = np.abs((d_dX2 - d_dX2Analytical)[domainX3Start:domainX3End, \
                                           domainX2Start:domainX2End, \
                                           domainX1Start:domainX1End])

errorX3 = np.abs((d_dX3 - d_dX3Analytical)[domainX3Start:domainX3End, \
                                           domainX2Start:domainX2End, \
                                           domainX1Start:domainX1End])

errorX1Left  = np.sum(errorX1[:, :, 0])  / (prim.N2Local*prim.N3Local)
errorX1Right = np.sum(errorX1[:, :, -1]) / (prim.N2Local*prim.N3Local)

errorX2Bottom = np.sum(errorX2[:, 0,  :]) / (prim.N1Local*prim.N3Local)
errorX2Top    = np.sum(errorX2[:, -1, :]) / (prim.N1Local*prim.N3Local)

errorX3Back  = np.sum(errorX3[0,  :, :])  / (prim.N1Local*prim.N2Local)
errorX3Front = np.sum(errorX3[-1, :, :])  / (prim.N1Local*prim.N2Local)

errorX1Bulk  = np.sum(errorX1[:, :, 1:-1]) / \
               ((prim.N1Local - 2) * prim.N3Local * prim.N1Local)

errorX2Bulk  = np.sum(errorX2[:, 1:-1, :]) / \
               (prim.N1Local * (prim.N2Local - 2) * prim.N3Local)

errorX3Bulk  = np.sum(errorX3[1:-1, :, :]) / \
               (prim.N1Local * prim.N2Local * (prim.N3Local - 2) )

print "----------------------------------"
print "rank = ", rank, " (i, j, k) = ", \
    "(", prim.iLocalStart, ",", prim.jLocalStart, ",", prim.kLocalStart, ")"

print "Error X1 left boundary   = ", errorX1Left
print "Error X1 right boundary  = ", errorX1Right
print "Error X1 bulk            = ", errorX1Bulk 
print ""
print "Error X2 bottom boundary = ", errorX2Bottom
print "Error X2 top boundary    = ", errorX2Top
print "Error X2 bulk            = ", errorX2Bulk
print ""
print "Error X3 back boundary   = ", errorX3Back
print "Error X3 front boundary  = ", errorX3Front
print "Error X3 bulk            = ", errorX3Bulk
print ""

# Error when communication fails is ~.6. Choose a value lower than that, but not
# too low that one would require a high grid resolution, since we are computing
# derivatives at O(dx^2)
errorTolerance = 1e-1

def test_communication_X1_left():
  assert np.abs(errorX1Left - errorX1Bulk) < errorTolerance

def test_communication_X1_right():
  assert np.abs(errorX1Right - errorX1Bulk) < errorTolerance

def test_communication_X2_top():
  assert np.abs(errorX2Top - errorX2Bulk) < errorTolerance

def test_communication_X2_bottom():
  assert np.abs(errorX2Bottom - errorX2Bulk) < errorTolerance

def test_communication_X3_front():
  assert np.abs(errorX3Front - errorX3Bulk) < errorTolerance

def test_communication_X3_back():
  assert np.abs(errorX3Back - errorX3Bulk) < errorTolerance

X1CoordsCheck = \
  (prim.iLocalStart + \
   np.arange(-numGhostX1, prim.N1Local + numGhostX1) + 0.5 \
  )*XCoords.dX1

X2CoordsCheck = \
  (prim.jLocalStart + \
   np.arange(-numGhostX2, prim.N2Local + numGhostX2) + 0.5 \
  )*XCoords.dX2

X3CoordsCheck = \
  (prim.kLocalStart + \
   np.arange(-numGhostX2, prim.N3Local + numGhostX3) + 0.5 \
  )*XCoords.dX3

X3CoordsCheck, X2CoordsCheck, X1CoordsCheck = \
    np.meshgrid(X3CoordsCheck, X2CoordsCheck, X1CoordsCheck, indexing='ij')

print "N1Local = ", prim.N1Local
print "N2Local = ", prim.N2Local
print "N3Local = ", prim.N3Local
print "X1Coords.shape      = ", X1Coords.shape
print "X1CoordsCheck.shape = ", X1CoordsCheck.shape
print "X2Coords.shape      = ", X2Coords.shape
print "X2CoordsCheck.shape = ", X2CoordsCheck.shape
print "X3Coords.shape      = ", X3Coords.shape
print "X3CoordsCheck.shape = ", X3CoordsCheck.shape

def test_X1Coords():
  assert np.sum(X1CoordsCheck - X1Coords) == 0

def test_X2Coords():
  assert np.sum(X2CoordsCheck - X2Coords) == 0

def test_X3Coords():
  assert np.sum(X3CoordsCheck - X3Coords) == 0
