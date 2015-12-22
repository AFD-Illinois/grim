#ifndef GRIM_GRID_H_
#define GRIM_GRID_H_

#include "../params.hpp"
#include <petsc.h>
#include <arrayfire.h>

using af::array;
using af::span;
using af::where;

class grid
{
  public:
    DM dm, dmGhost;
    Vec globalVec, localVec;

    int numVars, numGhost;

    DMBoundaryType boundaryLeft,  boundaryRight;
    DMBoundaryType boundaryTop,   boundaryBottom;
    DMBoundaryType boundaryFront, boundaryBack;

    array *vars;

    grid(int numVars, int numGhost);
    ~grid();

    void setVarsWithVec(Vec vec);
};

#endif /* GRIM_GRID_H_ */
