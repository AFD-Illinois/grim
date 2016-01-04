import gridPy
from gridPy cimport gridPy
from gridPythonHeaders cimport grid
from reconstructionPythonHeaders cimport reconstruct

def reconstructPy(gridPy prim, dir, gridPy primLeft, gridPy primRight):
#  primPtr = prim.getGridPtr()
#  cdef grid *primLeftPtr  = primLeft.getGridPtr()
#  cdef grid *primRightPtr = primRight.getGridPtr()
  reconstruct(prim.getGridPtr()[0], dir,
              primLeft.getGridPtr()[0], primRight.getGridPtr()[0]
             )

  return
