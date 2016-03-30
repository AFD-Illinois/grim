import mpi4py, petsc4py
from petsc4py import PETSc
import numpy as np
import pytest
import gridPy
import geometryPy
import boundaryPy
import timeStepperPy

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

Rin = 0.98*(1.+np.sqrt(1.-blackHoleSpin*blackHoleSpin));
Rout = 40.

X1Start = np.log(Rin); X1End = np.log(Rout)
X2Start = 1e-8; X2End = 1.-1e-8
X3Start = 0.; X3End = 2.*np.pi
periodicBoundariesX1 = False
periodicBoundariesX2 = False
periodicBoundariesX3 = False

boundaryLeft   = boundaryPy.OUTFLOW
boundaryRight  = boundaryPy.OUTFLOW
boundaryTop    = boundaryPy.OUTFLOW
boundaryBottom = boundaryPy.OUTFLOW
boundaryFront  = boundaryPy.PERIODIC
boundaryBack   = boundaryPy.PERIODIC


time = 0.
dt   = 0.01
numVars = 8
metric = geometryPy.MODIFIED_KERR_SCHILD
ts = timeStepperPy.timeStepperPy(N1, N2, N3,
                                 dim, numVars, numGhost,
                                 time, dt,
                                 boundaryLeft, boundaryRight,
                                 boundaryTop,  boundaryBottom,
                                 boundaryFront, boundaryBack,
                                 metric, blackHoleSpin, hSlope,
                                 X1Start, X1End,
                                 X2Start, X2End,
                                 X3Start, X3End
                                )
ts.computeDivB(ts.primOld)
print "divB.shape = ", ts.divB.shape
print "Div B = ", np.max(np.abs(ts.divB.getVars()[0, 0, numGhost:N2+numGhost,
  numGhost:N1+numGhost]))

import pylab as pl
for n in xrange(1000):
  print "Time step = ", n
  ts.timeStep()
  ts.computeDivB(ts.primHalfStep)
  print "Div B primHalfStep = ", \
    np.max(np.abs(ts.divB.getVars()[0, 0, numGhost:N2+numGhost, numGhost:N1+numGhost]))
  ts.computeDivB(ts.primOld)
  print "Div B primOld = ", \
    np.max(np.abs(ts.divB.getVars()[0, 0, numGhost:N2+numGhost, numGhost:N1+numGhost]))

  pl.contourf((np.abs(ts.divB.getVars()[0, 0, :, :])), 100)
  pl.colorbar()
  pl.savefig("divB_" + str(n) + ".png")
  pl.clf()
