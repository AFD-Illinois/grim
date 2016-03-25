#include "timestepper.hpp"

void timeStepper::fluxCT(grid &fluxX1,
                         grid &fluxX2,
                         grid &fluxX3,
                         int &numReads,
                         int &numWrites
                        )
{
  if (fluxX1.dim >= 2)
  {
    int numReadsEMF, numWritesEMF;
    computeEMF(fluxX1, fluxX2, fluxX3,
               *emfX1,  *emfX2,  *emfX3,
               numReadsEMF, numWritesEMF
              );

    fluxX1.vars[vars::B2] = 
      0.5*(  emfX3->vars[0]
           + shift(emfX3->vars[0], 0, -1, 0)
          );

    fluxX2.vars[vars::B1] =
      -0.5*(  emfX3->vars[0]
            + shift(emfX3->vars[0], -1, 0, 0)
           );
  
    if (fluxX1.dim == 3)
    {
      fluxX1.vars[vars::B3] =
        -0.5*(  emfX2->vars[0]
              + shift(emfX2->vars[0], 0, 0, -1)
             );

      fluxX2.vars[vars::B3] =
        0.5*(  emfX1->vars[0]
             + shift(emfX1->vars[0], 0, 0, -1)
            );
      
      fluxX3.vars[vars::B1] =
        0.5*(  emfX2->vars[0]
             + shift(emfX2->vars[0], -1, 0, 0)
            );

      fluxX3.vars[vars::B2] =
        -0.5*(  emfX1->vars[0]
              + shift(emfX1->vars[0], 0, -1, 0)
             );
    }
  }




}

void timeStepper::computeEMF(const grid &fluxX1,
                             const grid &fluxX2,
                             const grid &fluxX3,
                             grid &emfX1,
                             grid &emfX2,
                             grid &emfX3,
                             int &numReadsEMF,
                             int &numWritesEMF
                            )
{
  if (fluxX1.dim >= 2)
  {
    emfX3.vars[0] = 
      0.25*(  fluxX1.vars[vars::B2] 
            + shift(fluxX1.vars[vars::B2], 0, 1, 0)
            - fluxX2.vars[vars::B1]
            - shift(fluxX2.vars[vars::B1], 1, 0, 0)
           );


    if (fluxX1.dim == 3)
    {
      emfX1.vars[0] =
        0.25*(  fluxX2.vars[vars::B3]
              + shift(fluxX2.vars[vars::B3], 0, 0, 1)
              - fluxX3.vars[vars::B2]
              - shift(fluxX3.vars[vars::B2], 0, 1, 0)
             );

      emfX2.vars[0] =
        0.25*(  fluxX3.vars[vars::B1]
              + shift(fluxX3.vars[vars::B1], 1, 0, 0)
              - fluxX1.vars[vars::B3]
              - shift(fluxX1.vars[vars::B3], 0, 0, 1)
             );
    }
  }
}

void timeStepper::computeDivB(const grid &prim,
                              grid &divB,
                              int &numReads,
                              int &numWrites
                             )
{

  array B1 = prim.vars[vars::B1];
  array B2 = prim.vars[vars::B2];
  array B3 = prim.vars[vars::B3];
  array g  = geomCenter->g;

  divB.vars[0] = 
      (  shift(g*B1, 0, -1, 0) + shift(g*B1, 0, -1, -1)
       - g*B1 - shift(g*B1, 0, 0, 1)
      )/(2.*XCoords->dX1)
   + 
      (  shift(g*B2, 0, 0, -1) + shift(g*B2, 0, -1, -1)
       - g*B2 - shift(g*B2, 0, -1, 0)
      );
}
