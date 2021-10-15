#include "timestepper.hpp"

void timeStepper::computeDivOfFluxes(const grid &primFlux,
                                     int &numReads,
                                     int &numWrites
                                    )
{
  int numReadsReconstruction, numWritesReconstruction;
  int numReadsRiemann, numWritesRiemann;
  int numReadsCT, numWritesCT;

  switch (primFlux.dim)
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

      applyProblemSpecificFluxFilter(numReads,numWrites);

      for (int var=0; var < primFlux.numVars; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(XCoords->dX1);

        array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);

        divFluxes->vars[var] = dFluxX1_dX1;
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

      fluxCT(numReadsCT, numWritesCT);
      numReads  += numReadsCT;
      numWrites += numWritesCT;

      applyProblemSpecificFluxFilter(numReads,numWrites);

      for (int var=0; var < primFlux.numVars; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(XCoords->dX1);
        array filterX2 = array(1, 3, 1, 1, filter1D)/(XCoords->dX2);

        array dFluxX1_dX1 = convolve(fluxesX1->vars[var], filterX1);
        array dFluxX2_dX2 = convolve(fluxesX2->vars[var], filterX2);
        
        divFluxes->vars[var] = dFluxX1_dX1 + dFluxX2_dX2;
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

      fluxCT(numReadsCT, numWritesCT);
      numReads  += numReadsCT;
      numWrites += numWritesCT;

      applyProblemSpecificFluxFilter(numReads,numWrites);

      for (int var=0; var < primFlux.numVars; var++)
      {
        array dFluxX1_dX1 = (  af::shift(fluxesX1->vars[var], -1,  0,  0) 
                             - fluxesX1->vars[var]
                            )/XCoords->dX1;

        array dFluxX2_dX2 = (  af::shift(fluxesX2->vars[var],  0, -1,  0)
                             - fluxesX2->vars[var]
                            )/XCoords->dX2;

        array dFluxX3_dX3 = (  af::shift(fluxesX3->vars[var],  0,  0, -1)
                             - fluxesX3->vars[var]
                            )/XCoords->dX3;

        divFluxes->vars[var] = dFluxX1_dX1 + dFluxX2_dX2 + dFluxX3_dX3;
      }

      break;
  }
}

