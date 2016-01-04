from gridPy cimport gridPy # cimport stands for importing C types. Needed here 
                           # to set the definitions for the arguments of
                           # reconstructPy
from reconstructionPythonHeaders cimport reconstruct # Import the C function

def reconstructPy(gridPy prim, dir, gridPy primLeft, gridPy primRight):
  reconstruct(prim.getGridPtr()[0], dir,
              primLeft.getGridPtr()[0], primRight.getGridPtr()[0]
             )

  return
