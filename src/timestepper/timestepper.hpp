#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

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
  grid *primGuessLineSearchTrial;
  grid *primGuessPlusEps;
  grid *residual;
  grid *residualPlusEps;

  array residualSoA;
  array jacobianSoA;
  array deltaPrimAoS;
  array stepLength;

  double *AHostPtr, *bHostPtr, *xHostPtr;

  void solve(grid &primGuess);
  void computeResidual(const grid &prim, grid &residual,
		                   const bool ComputeExplicitTerms,
                       int &numReads,
                       int &numWrites
                      );
  void batchLinearSolve(const array &A, const array &b, array &x);

  public:
    double dt, time;
    int N1, N2, N3;
    int numVars;

    grid *XCoords;
    grid *prim, *primHalfStep, *primOld;
    grid *cons, *consOld;
    grid *sourcesExplicit;
    grid *sourcesImplicit;
    grid *sourcesImplicitOld;
    grid *sourcesTimeDer;
    grid *primLeft, *primRight;
    grid *fluxesX1, *fluxesX2, *fluxesX3;
    grid *divFluxes;

    geometry *geomLeft, *geomRight;
    geometry *geomBottom, *geomTop;
    geometry *geomCenter;

    fluidElement *elem, *elemOld, *elemHalfStep;

    riemannSolver *riemann;

    void computeDivOfFluxes(const grid &prim,
                            int &numReads, int &numWrites
                           );

    int currentStep;

    timeStepper(const int N1, const int N2, const int N3,
                const int numGhost, const int dim, 
                const int numVars, 
                const double time,
                const double dt,
                DMBoundaryType boundaryLeft,  DMBoundaryType boundaryRight,
                DMBoundaryType boundaryTop,   DMBoundaryType boundaryBottom,
                DMBoundaryType boundaryFront, DMBoundaryType boundaryBack
               );
    ~timeStepper();

    void timeStep(int &numReads, int &numWrites);

    /* Function definitions in the problem folder */
    void initialConditions(const array xCoords[3], array prim[vars::dof]);
    void halfStepDiagnostics(const array xCoords[3], array prim[vars::dof]);
    void fullStepDiagnostics(const array xCoords[3], array prim[vars::dof]);
    void setProblemSpecificBCs(const array xCoords[3], array prim[vars::dof]);
};

#endif /* GRIM_TIMESTEPPER_H_ */
