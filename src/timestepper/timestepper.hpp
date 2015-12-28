#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../physics/physics.hpp"
#include "../geometry/geometry.hpp"
#include "../boundary/boundary.hpp"
#include "../nonlinearsolver/nonlinearsolver.hpp"

namespace timeStepperSwitches
{
  enum
  {
    HALF_STEP, FULL_STEP
  };
};

class timeStepper
{
  public:
    geometry *geom;

    grid *prim, *primHalfStep, *primOld;
    grid *cons, *consOld;
    grid *sources, *sourcesHalfStep, *sourcesOld;
    grid *fluxesX1, *fluxesX2, *fluxesX3;
    grid *divFluxes;

    fluidElement *elem, *elemOld, *elemHalfStep;

    riemannSolver *riemann;
    nonLinearSolver *nonLinSolver;

    void computeDivOfFluxes(const grid &primGhosted);

    int currentStep;

    timeStepper();
    ~timeStepper();

    void timeStep(double dt);
};

void computeResidual(const grid &prim, grid &residual, void *dataPtr);

#endif /* GRIM_TIMESTEPPER_H_ */
