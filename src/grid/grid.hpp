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

    int numGhost, numVars, dim;

    DMBoundaryType boundaryLeft,  boundaryRight;
    DMBoundaryType boundaryTop,   boundaryBottom;
    DMBoundaryType boundaryFront, boundaryBack;

    int iLocalStart, jLocalStart, kLocalStart;
    int iLocalEnd,   jLocalEnd,   kLocalEnd;
    int N1,          N2,          N3;
    int N1Local,     N2Local,     N3Local;
    int N1Total,     N2Total,     N3Total;
    int numGhostX1,  numGhostX2,  numGhostX3;

    af::seq *domainX1, *domainX2, *domainX3;

    double dX1, dX2, dX3;

    array *vars, varsSoA;
    array indices[3];
    double *hostPtr;
    bool hasHostPtrBeenAllocated;

    grid(int N1, int N2, int N3, 
         int numGhost, int dim, int numVars,
         DMBoundaryType boundaryLeft,
         DMBoundaryType boundaryRight,
         DMBoundaryType boundaryTop,
         DMBoundaryType boundaryBottom,
         DMBoundaryType boundaryFront,
         DMBoundaryType boundaryBack
        );
    ~grid();

    void communicate();
    void copyVarsToHostPtr();
    void copyHostPtrToVars(const double *hostPtr);
};

#endif /* GRIM_GRID_H_ */
