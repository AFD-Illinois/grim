#ifndef GRIM_GRID_H_
#define GRIM_GRID_H_

#include "../params.hpp"
#include <petsc.h>
#include <arrayfire.h>

using af::array;
using af::span;

class grid
{
  void copyLocalVecToVars();

  public:
    DM dm;
    Vec globalVec, localVec;

    int numVars, numGhost;

    DMBoundaryType boundaryLeft,  boundaryRight;
    DMBoundaryType boundaryTop,   boundaryBottom;
    DMBoundaryType boundaryFront, boundaryBack;

    int iLocalStart, jLocalStart, kLocalStart;
    int iLocalEnd,   jLocalEnd,   kLocalEnd;
    int N1Local,     N2Local,     N3Local;
    int N1Total,     N2Total,     N3Total;
    int numGhostX1,  numGhostX2,  numGhostX3;

    af::seq *domainX1, *domainX2, *domainX3;

    double dX1, dX2, dX3;

    array *vars, varsSoA;

    grid(int numVars);
    ~grid();

    void communicate();
};

#endif /* GRIM_GRID_H_ */
