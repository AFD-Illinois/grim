from gridPy cimport gridPy # cimport stands for importing C types. Needed here 
                           # to set the definitions for the arguments of
                           # reconstructPy
# Import the C++ functions
from reconstructionPythonHeaders cimport reconstruct
from reconstructionPythonHeaders cimport reconstructMM
from reconstructionPythonHeaders cimport reconstructWENO5

def reconstructPy(gridPy prim, dir, gridPy primLeft, gridPy primRight):
  reconstruct(prim.getGridPtr()[0], dir,
              primLeft.getGridPtr()[0], primRight.getGridPtr()[0]
             )
  return

def reconstructMinModPy(gridPy prim, dir, gridPy primLeft, gridPy primRight):
  reconstructMM(prim.getGridPtr()[0], dir,
                primLeft.getGridPtr()[0], primRight.getGridPtr()[0]
               )
  return

def reconstructWENO5Py(gridPy prim, dir, gridPy primLeft, gridPy primRight):
  reconstructWENO5(prim.getGridPtr()[0], dir,
                   primLeft.getGridPtr()[0], primRight.getGridPtr()[0]
                  )
  return
