import mpi4py, petsc4py
import numpy as np
import pytest
import gridPy
import boundaryPy

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
periodicBoundariesX1 = False
periodicBoundariesX2 = False
periodicBoundariesX3 = False

boundaryLeft   = boundaryPy.OUTFLOW
boundaryRight  = boundaryPy.OUTFLOW
boundaryTop    = boundaryPy.OUTFLOW
boundaryBottom = boundaryPy.OUTFLOW
boundaryFront  = boundaryPy.OUTFLOW
boundaryBack   = boundaryPy.OUTFLOW

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

varsNumpy = np.random.rand(prim.shape[0], 
                           prim.shape[1],
                           prim.shape[2],
                           prim.shape[3]
                          )

prim.setVars(varsNumpy)

def test_outflow_X1Left():
  boundaryLeft   = boundaryPy.OUTFLOW
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.iLocalStart==0):
    for var in xrange(prim.numVars):
      for i in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, :, i], 
                                varsNumpy[var, :, :, numGhost]
                               )
def test_outflow_X1Right():
  boundaryRight   = boundaryPy.OUTFLOW
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.iLocalEnd==N1):
    for var in xrange(prim.numVars):
      for i in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, :, prim.N1Local+numGhost+i], 
                                varsNumpy[var, :, :, prim.N1Local+numGhost-1]
                               )

def test_mirror_X1Left():
  boundaryLeft   = boundaryPy.MIRROR
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.iLocalStart==0):
    for var in xrange(prim.numVars):
      for i in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, :, i], 
                                varsNumpy[var, :, :, 2*numGhost-i-1]
                               )

def test_mirror_X1Right():
  boundaryRight   = boundaryPy.MIRROR
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.iLocalEnd==N1):
    for var in xrange(prim.numVars):
      for i in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, :, prim.N1Local+numGhost+i], 
                                varsNumpy[var, :, :, prim.N1Local+numGhost-1-i]
                               )

def test_outflow_X2Bottom():
  boundaryBottom = boundaryPy.OUTFLOW
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.jLocalStart==0):
    for var in xrange(prim.numVars):
      for j in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, j, :], 
                                varsNumpy[var, :, numGhost, :]
                               )
def test_outflow_X2Top():
  boundaryTop = boundaryPy.OUTFLOW
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.jLocalEnd==N2):
    for var in xrange(prim.numVars):
      for j in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, prim.N2Local+numGhost+j, :],
                                varsNumpy[var, :, prim.N2Local+numGhost-1, :]
                               )

def test_mirror_X2Bottom():
  boundaryBottom = boundaryPy.MIRROR
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.jLocalStart==0):
    for var in xrange(prim.numVars):
      for j in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, j, :], 
                                varsNumpy[var, :, 2*numGhost-j-1, :]
                               )

def test_mirror_X2Top():
  boundaryTop = boundaryPy.MIRROR
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.jLocalEnd==N2):
    for var in xrange(prim.numVars):
      for j in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, :, prim.N2Local+numGhost+j, :],
                                varsNumpy[var, :, prim.N2Local+numGhost-1-j, :]
                               )

def test_outflow_X3Back():
  boundaryBack = boundaryPy.OUTFLOW
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.kLocalStart==0):
    for var in xrange(prim.numVars):
      for k in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, k, :, :], 
                                varsNumpy[var, numGhost, :, :]
                               )
def test_outflow_X3Front():
  boundaryFront = boundaryPy.OUTFLOW
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.kLocalEnd==N3):
    for var in xrange(prim.numVars):
      for k in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, prim.N3Local+numGhost+k, :, :],
                                varsNumpy[var, prim.N3Local+numGhost-1, :, :]
                               )

def test_mirror_X3Back():
  boundaryBack = boundaryPy.MIRROR
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.kLocalStart==0):
    for var in xrange(prim.numVars):
      for k in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, k, :, :], 
                                varsNumpy[var, 2*numGhost-k-1, :, :]
                               )

def test_mirror_X3Front():
  boundaryFront = boundaryPy.MIRROR
  boundaryPy.applyBoundaryConditionsPy(boundaryLeft,  boundaryRight,
                                       boundaryTop,   boundaryBottom,
                                       boundaryFront, boundaryBack,
                                       prim
                                      )
  varsNumpy = prim.getVars()

  if (prim.kLocalEnd==N3):
    for var in xrange(prim.numVars):
      for k in xrange(numGhost):
        np.testing.assert_equal(varsNumpy[var, prim.N3Local+numGhost+k, :, :],
                                varsNumpy[var, prim.N3Local+numGhost-1-k, :, :]
                               )
