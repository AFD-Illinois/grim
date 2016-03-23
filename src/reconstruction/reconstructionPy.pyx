from gridPy cimport gridPy # cimport stands for importing C types. Needed here 
                           # to set the definitions for the arguments of
                           # reconstructPy
# Import the C++ functions
from reconstructionHeaders cimport reconstruct
from reconstructionHeaders cimport reconstructMM
from reconstructionHeaders cimport reconstructWENO5

def reconstructPy(gridPy prim,
                  int dir, 
                  gridPy primLeft,
                  gridPy primRight
                 ):
  cdef int numReads = 0
  cdef int numWrites = 0
  reconstruct(prim.getGridPtr()[0], 
              dir,
              primLeft.getGridPtr()[0],
              primRight.getGridPtr()[0],
              numReads, numWrites
             )
  return numReads, numWrites

def reconstructMinModPy(gridPy prim, 
                        int dir,
                        gridPy primLeft,
                        gridPy primRight
                       ):
  cdef int numReads  = 0
  cdef int numWrites = 0
  reconstructMM(prim.getGridPtr()[0],
                dir,
                primLeft.getGridPtr()[0],
                primRight.getGridPtr()[0],
                numReads, numWrites
               )
  return numReads, numWrites

def reconstructWENO5Py(gridPy prim, dir, gridPy primLeft, gridPy primRight):
  cdef int numReads  = 0
  cdef int numWrites = 0
  reconstructWENO5(prim.getGridPtr()[0],
                   dir,
                   primLeft.getGridPtr()[0],
                   primRight.getGridPtr()[0],
                   numReads, numWrites
                  )
  return numReads, numWrites
