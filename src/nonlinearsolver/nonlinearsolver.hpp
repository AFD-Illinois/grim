#ifndef GRIM_NONLINEARSOLVER_H_
#define GRIM_NONLINEARSOLVER_H_

#include "../params.hpp"
#include "../grid/grid.hpp"

class nonLinearSolver
{
  grid *primGuessLineSearchTrial;
  grid *primGuessPlusEps;
  grid *residual;
  grid *residualPlusEps;

  array residualSoA;
  array jacobianSoA;
  array bSoA;
  array deltaPrimSoA;
  array stepLength;

  void (*computeResidual)(const grid &primGuess, grid &residual);

  public:

    nonLinearSolver(void (*computeResidual)(const grid &primGuess, 
                                            grid &residual
                                           )
                   );
    ~nonLinearSolver();

    void solve(grid &primGuess);
};

#endif /* GRIM_NONLINEARSOLVER_H_ */
