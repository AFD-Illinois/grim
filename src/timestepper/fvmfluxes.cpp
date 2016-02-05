#include "timestepper.hpp"

void timeStepper::computeDivOfFluxes(const grid &primFlux,
                                     int &numReads,
                                     int &numWrites
                                    )
{
  int numReadsReconstruction, numWritesReconstruction;
  int numReadsRiemann, numWritesRiemann;

  switch (params::dim)
  {
    case 1:
    /* Reconstruction gives, at a point of index i:
     * primLeft : right-biased stencil reconstructs on face i-/1.2
     * primRight: left-biased stencil reconstructs on face i+1/2 */
      reconstruction::reconstruct(primFlux, directions::X1,
                                  *primLeft, *primRight,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  = numReadsReconstruction;
      numWrites = numWritesReconstruction;

      riemann->solve(*primLeft, *primRight,
                     *geomLeft, *geomRight,
                     directions::X1, *fluxesX1,
                     numReadsRiemann, numWritesRiemann
                    );
      numReads  += numReadsRiemann;
      numWrites += numWritesRiemann;

      for (int var=0; var < primFlux.numVars; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxesX1->dX1);

        array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);

        divFluxes->vars[var] = dFluxX1_dX1;
        divFluxes->vars[var].eval();
      }

      break;

    case 2:

      /* directions:: X1 */
      reconstruction::reconstruct(primFlux, directions::X1,
                                  *primLeft, *primRight,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  = numReadsReconstruction;
      numWrites = numWritesReconstruction;

      riemann->solve(*primLeft, *primRight,
                     *geomLeft, *geomRight,
                     directions::X1, *fluxesX1,
                     numReadsRiemann, numWritesRiemann
                    );
      numReads  += numReadsRiemann;
      numWrites += numWritesRiemann;

      /* directions:: X2 */
      reconstruction::reconstruct(primFlux, directions::X2,
                                  *primLeft, *primRight,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  += numReadsReconstruction;
      numWrites += numWritesReconstruction;

      riemann->solve(*primLeft,   *primRight,
                     *geomBottom, *geomTop,
                     directions::X2, *fluxesX2,
                     numReadsRiemann, numWritesRiemann
                    );
      numReads  += numReadsRiemann;
      numWrites += numWritesRiemann;

      for (int var=0; var < primFlux.numVars; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxesX1->dX1);
        array filterX2 = array(1, 3, 1, 1, filter1D)/(fluxesX2->dX2);

        array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);
        array dFluxX2_dX2 = convolve(fluxesX2->vars[var], filterX2);

        divFluxes->vars[var] = dFluxX1_dX1 + dFluxX2_dX2;
        divFluxes->vars[var].eval();
      }

      break;

    case 3:
      /* directions:: X1 */
      reconstruction::reconstruct(primFlux, directions::X1,
                                  *primLeft, *primRight,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  = numReadsReconstruction;
      numWrites = numWritesReconstruction;

      riemann->solve(*primLeft, *primRight,
                     *geomLeft, *geomRight,
                     directions::X1, *fluxesX1,
                     numReadsRiemann, numWritesRiemann
                    );
      numReads  += numReadsRiemann;
      numWrites += numWritesRiemann;

      /* directions:: X2 */
      reconstruction::reconstruct(primFlux, directions::X2,
                                  *primLeft, *primRight,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  += numReadsReconstruction;
      numWrites += numWritesReconstruction;

      riemann->solve(*primLeft,   *primRight,
                     *geomBottom, *geomTop,
                     directions::X2, *fluxesX2,
                     numReadsRiemann, numWritesRiemann
                    );
      numReads  += numReadsRiemann;
      numWrites += numWritesRiemann;

      /* directions:: X3 */
      reconstruction::reconstruct(primFlux, directions::X3,
                                  *primLeft, *primRight,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  += numReadsReconstruction;
      numWrites += numWritesReconstruction;

      riemann->solve(*primLeft,   *primRight,
                     *geomCenter, *geomCenter,
                     directions::X3, *fluxesX3,
                     numReadsRiemann, numWritesRiemann
                    );
      numReads  += numReadsRiemann;
      numWrites += numWritesRiemann;

      for (int var=0; var < primFlux.numVars; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxesX1->dX1);
        array filterX2 = array(1, 3, 1, 1, filter1D)/(fluxesX2->dX2);
        array filterX3 = array(1, 1, 3, 1, filter1D)/(fluxesX3->dX3);

        array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);
        array dFluxX2_dX2 = convolve(fluxesX2->vars[var], filterX2);
        array dFluxX3_dX3 = convolve(fluxesX3->vars[var], filterX3);

        divFluxes->vars[var] = dFluxX1_dX1 + dFluxX2_dX2 + dFluxX3_dX3;
        divFluxes->vars[var].eval();
      }

      break;
  }
}

