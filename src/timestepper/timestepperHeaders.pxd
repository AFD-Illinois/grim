from gridHeaders cimport grid, coordinatesGrid

cdef extern from "timestepper.hpp":
  cdef enum:
    STAGE_HALF_STEP "timeStepperSwitches::HALF_STEP"
    STAGE_FULL_STEP "timeStepperSwitches::FULL_STEP"

  cdef cppclass timeStepper:
    timeStepper(const int N1, 
                const int N2,
                const int N3,
                const int dim,
                const int numVars, 
                const int numGhost,
                const double time,
                const double dt,
                const int boundaryLeft,  const int boundaryRight,
                const int boundaryTop,   const int boundaryBottom,
                const int boundaryFront, const int boundaryBack,
                const int metric,
                const double blackHoleSpin,
                const double hSlope,
                const double X1Start, const double X1End,
                const double X2Start, const double X2End,
                const double X3Start, const double X3End
               )
    coordinatesGrid *XCoords
    grid *prim
    grid *primOld
    grid *primHalfStep
    grid *fluxesX1
    grid *fluxesX2
    grid *fluxesX3
    grid *cons
    grid *consOld
    grid *sourcesExplicit
    grid *sourcesImplicit
    grid *sourcesImplicitOld
    grid *sourcesTimeDer

    void timeStep(int &numReads, int &numWrites)

    void fluxCT(int &numReads, int &numWrites)
    void computeEMF(int &numReadsEMF, int &numWritesEMF)
    void computeDivB(const grid &prim,
                     grid &divB,
                     int &numReads,
                     int &numWrites
                    )
