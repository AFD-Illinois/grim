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
  array deltaPrimAoS;
  array stepLength;

  void *dataPtr;

  void (*computeResidual)(const grid &primGuess, grid &residual, void *dataPtr);

  public:

    nonLinearSolver(void (*computeResidual)(const grid &primGuess, 
                                            grid &residual,
                                            void *dataPtr
                                           )
                   );
    ~nonLinearSolver();

    void solve(grid &primGuess);
};

#endif /* GRIM_NONLINEARSOLVER_H_ */
