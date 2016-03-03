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
         
    int boundaryLeft,  boundaryRight;
    int boundaryTop,   boundaryBottom;
    int boundaryFront, boundaryBack;
    int periodicBoundariesX1;
    int periodicBoundariesX2;
    int periodicBoundariesX3;

    DMBoundaryType DMBoundaryLeft,  DMBoundaryRight;
    DMBoundaryType DMBoundaryTop,   DMBoundaryBottom;
    DMBoundaryType DMBoundaryFront, DMBoundaryBack;

    int iLocalStart, jLocalStart, kLocalStart;
    int iLocalEnd,   jLocalEnd,   kLocalEnd;
    int N1,          N2,          N3;
    int N1Local,     N2Local,     N3Local;
    int N1Total,     N2Total,     N3Total;
    int numGhostX1,  numGhostX2,  numGhostX3;

    af::seq *domainX1, *domainX2, *domainX3;

    array *vars, varsSoA;
    array indices[3];
    double *hostPtr;
    bool hasHostPtrBeenAllocated;

    grid(int N1, int N2, int N3, 
         int dim, 
         int numVars,
         int numGhost,
         int periodicBoundariesX1,
         int periodicBoundariesX2,
         int periodicBoundariesX3
        );
    ~grid();

    void communicate();
    void copyVarsToHostPtr();
    void copyHostPtrToVars(const double *hostPtr);
};

class coordinatesGrid : public grid
{
  public:
    double X1Start, X1End;
    double X2Start, X2End;
    double X3Start, X3End;

    double dX1, dX2, dX3;

    coordinatesGrid(int N1, int N2, int N3,
                    int dim,
                    int numGhost,
                    double X1Start, double X1End,
                    double X2Start, double X2End,
                    double X3Start, double X3End
                   );
    ~coordinatesGrid();

    void setXCoords(const int location);

};

#endif /* GRIM_GRID_H_ */
