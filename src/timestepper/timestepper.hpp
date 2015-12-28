#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../physics/physics.hpp"
#include "../geometry/geometry.hpp"
#include "../boundary/boundary.hpp"

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
  array bSoA;
  array deltaPrimAoS;
  array stepLength;

  void solve(grid &primGuess);
  void computeResidual(const grid &prim, grid &residual);

  public:
    geometry *geom;

    grid *prim, *primHalfStep, *primOld;
    grid *cons, *consOld;
    grid *sources, *sourcesHalfStep, *sourcesOld;
    grid *fluxesX1, *fluxesX2, *fluxesX3;
    grid *divFluxes;

    fluidElement *elem, *elemOld, *elemHalfStep;

    riemannSolver *riemann;

    void computeDivOfFluxes(const grid &primGhosted);

    int currentStep;

    timeStepper();
    ~timeStepper();

    void timeStep(double dt);
};

#endif /* GRIM_TIMESTEPPER_H_ */
