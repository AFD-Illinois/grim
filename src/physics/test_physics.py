import mpi4py, petsc4py
from petsc4py import PETSc
import numpy as np
import pytest
import gridPy
import geometryPy
import physicsPy

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
#X1Start = 0.; X1End = 1.
#X2Start = 0.; X2End = 1.
#X3Start = 0.; X3End = 1.
Rin = 0.98*(1.+np.sqrt(1.-blackHoleSpin*blackHoleSpin));
Rout = 40.

X1Start = np.log(Rin); X1End = np.log(Rout)
X2Start = 1e-8; X2End = 1.-1e-8
X3Start = 0.; X3End = 2.*np.pi
periodicBoundariesX1 = False
periodicBoundariesX2 = False
periodicBoundariesX3 = False

numVars=10
prim = gridPy.gridPy(N1, N2, N3, 
                     dim, numVars, numGhost,
                     periodicBoundariesX1,
                     periodicBoundariesX2,
                     periodicBoundariesX3
                    )
fluxesX1 = gridPy.gridPy(N1, N2, N3, 
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

geom = geometryPy.geometryPy(geometryPy.MODIFIED_KERR_SCHILD,
                             blackHoleSpin, hSlope,
                             XCoords
                            )
elem = physicsPy.fluidElementPy(prim, geom)
numReads, numWrites = elem.computeFluxes(gridPy.X1, fluxesX1)
print "numReads = ", numReads, " numWrites = ", numWrites
