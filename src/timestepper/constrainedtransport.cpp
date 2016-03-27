#include "timestepper.hpp"

void timeStepper::fluxCT(int &numReads,
                         int &numWrites
                        )
{
  if (fluxesX1->dim >= 2)
  {
    int numReadsEMF, numWritesEMF;
    computeEMF(numReadsEMF, numWritesEMF);

    fluxesX1->vars[vars::B2] = 
      0.5*(  emfX3->vars[0]
           + shift(emfX3->vars[0], 0, -1, 0)
          );

    fluxesX2->vars[vars::B1] =
      -0.5*(  emfX3->vars[0]
            + shift(emfX3->vars[0], -1, 0, 0)
           );
  
    if (fluxesX1->dim == 3)
    {
      fluxesX1->vars[vars::B3] =
        -0.5*(  emfX2->vars[0]
              + shift(emfX2->vars[0], 0, 0, -1)
             );

      fluxesX2->vars[vars::B3] =
        0.5*(  emfX1->vars[0]
             + shift(emfX1->vars[0], 0, 0, -1)
            );
      
      fluxesX3->vars[vars::B1] =
        0.5*(  emfX2->vars[0]
             + shift(emfX2->vars[0], -1, 0, 0)
            );

      fluxesX3->vars[vars::B2] =
        -0.5*(  emfX1->vars[0]
              + shift(emfX1->vars[0], 0, -1, 0)
             );
    }
  }

}

void timeStepper::computeEMF(int &numReadsEMF,
                             int &numWritesEMF
                            )
{
  if (fluxesX1->dim >= 2)
  {
    emfX3->vars[0] = 
      0.25*(  fluxesX1->vars[vars::B2] 
            + shift(fluxesX1->vars[vars::B2], 0, 1, 0)
            - fluxesX2->vars[vars::B1]
            - shift(fluxesX2->vars[vars::B1], 1, 0, 0)
           );


    if (fluxesX1->dim == 3)
    {
      emfX1->vars[0] =
        0.25*(  fluxesX2->vars[vars::B3]
              + shift(fluxesX2->vars[vars::B3], 0, 0, 1)
              - fluxesX3->vars[vars::B2]
              - shift(fluxesX3->vars[vars::B2], 0, 1, 0)
             );

      emfX2->vars[0] =
        0.25*(  fluxesX3->vars[vars::B1]
              + shift(fluxesX3->vars[vars::B1], 1, 0, 0)
              - fluxesX1->vars[vars::B3]
              - shift(fluxesX1->vars[vars::B3], 0, 0, 1)
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
