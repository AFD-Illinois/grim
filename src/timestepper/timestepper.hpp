#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../physics/physics.hpp"
#include "../geometry/geometry.hpp"
#include "../boundary/boundary.hpp"

class timeStepper
{
  public:
    SNES snes;
    grid *prim, *primOld;
    grid *cons, *consOld;
    grid *sourcesOld;
    grid *residual;

    grid *fluxesX1, *fluxesX2, *fluxesX3;

    fluidElement *elem, *elemOld;
    riemannSolver *riemann;

    timeStepper(const geometry &geom);
    ~timeStepper();

    void timeStep();
};

void consToPrim(grid &cons, geometry &geom, grid &primGuess);

#endif /* GRIM_TIMESTEPPER_H_ */
