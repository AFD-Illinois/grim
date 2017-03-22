#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include <sys/stat.h>
#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../physics/physics.hpp"
#include "../geometry/geometry.hpp"
#include "../boundary/boundary.hpp"
#include "mkl.h"

namespace timeStepperSwitches
{
  enum
  {
    HALF_STEP, FULL_STEP
  };
};

class timeStepper
{
  int world_rank, world_size;

  grid *primGuessLineSearchTrial;
  grid *primGuessPlusEps;
  grid *primIC;
  grid *residual;
  grid *residualPlusEps;
  grid *dump;

  array residualSoA;
  array jacobianSoA;
  array deltaPrimAoS;
  array stepLength;

  double *AHostPtr, *bHostPtr;

  void solve(grid &primGuess);
  void computeResidual(const grid &prim, grid &residual,
                       int &numReads,
                       int &numWrites
                      );
  void batchLinearSolve(const array &A, const array &b, array &x);
  double linearSolverTime;
  double lineSearchTime;
  double jacobianAssemblyTime;

  af::seq domainX1, domainX2, domainX3;
  array residualMask;

  double memoryBandwidth(const double numReads,
                         const double numWrites,
                         const double numEvals,
                         const double timeElapsed
                        );

  double bandwidthTest(const int numEvals);

  public:
    double dt, time;
    int N1, N2, N3, numGhost, dim;
    int numVars;

    int boundaryLeft, boundaryRight;
    int boundaryTop,  boundaryBottom;
    int boundaryFront, boundaryBack;

    coordinatesGrid *XCoords;
    grid *prim, *primHalfStep, *primOld;
    grid *cons, *consOld;
    grid *sourcesExplicit;
    grid *sourcesImplicit;
    grid *sourcesImplicitOld;
    grid *sourcesTimeDer;
    grid *primLeft, *primRight;
    grid *fluxesX1, *fluxesX2, *fluxesX3;
    grid *emfX1, *emfX2, *emfX3;
    grid *divFluxes;
    grid *divB;

    geometry *geomLeft,   *geomRight;
    geometry *geomBottom, *geomTop;
    geometry *geomCenter;

    fluidElement *elem, *elemOld, *elemHalfStep;

    riemannSolver *riemann;

    void computeDivOfFluxes(const grid &prim,
                            int &numReads, int &numWrites
                           );

    int currentStep;

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
               );
    ~timeStepper();

    void timeStep(int &numReads, int &numWrites);

    void fluxCT(int &numReads, int &numWrites);
    void computeEMF(int &numReadsEMF, int &numWritesEMF);
    void computeDivB(const grid &prim,
                     int &numReads,
                     int &numWrites
                    );

    double computeDt(int &numReads, int &numWrites);

    /* Function definitions in the problem folder */
    void initialConditions(int &numReads, int &numWrites);
    void halfStepDiagnostics(int &numReads, int &numWrites);
    void fullStepDiagnostics(int &numReads, int &numWrites);
    void setProblemSpecificBCs(int &numReads, int &numWrites);
    void applyProblemSpecificFluxFilter(int &numReads, int &numWrites);
    int CheckWallClockTermination();
};

#endif /* GRIM_TIMESTEPPER_H_ */
