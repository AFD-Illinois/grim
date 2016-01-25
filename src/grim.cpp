#include "grim.hpp"
#include "params.cpp"

void tiledComputation(const grid &prim, 
                      const geometry &geom, 
                      grid &cons,
                      fluidElement &elemTile,
                      geometry &geomTile,
                      array primTile[vars::dof],
                      array consTile[vars::dof]
                     );

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscPrintf(PETSC_COMM_WORLD, 
              "#### System info ####\n"
             );
  if (rank==0)
  {
    af::info();
  };

  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    grid indices(params::N1, params::N2, params::N3,
                 params::numGhost, params::dim, 3,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                );
    grid XCoords(params::N1, params::N2, params::N3,
                 params::numGhost, params::dim, 3,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                );

    indices.vars[directions::X1]
      = af::range(indices.N1Total, /* number of total zones in X1 */
                  indices.N2Total, /* number of total zones in X2 */
                  indices.N3Total, /* number of total zones in X3 */
                  1,                /* number of variables */
                  directions::X1 ,  /* Vary indices in X1 direction */
                  f64               /* Double precision */
                 ) - indices.numGhostX1;

    indices.vars[directions::X2]
      = af::range(indices.N1Total, /* number of total zones in X1 */
                  indices.N2Total, /* number of total zones in X2 */
                  indices.N3Total, /* number of total zones in X3 */
                  1,                /* number of variables */
                  directions::X2 ,  /* Vary indices in X1 direction */
                  f64               /* Double precision */
                 ) - indices.numGhostX2;

    indices.vars[directions::X3]
      = af::range(indices.N1Total, /* number of total zones in X1 */
                  indices.N2Total, /* number of total zones in X2 */
                  indices.N3Total, /* number of total zones in X3 */
                  1,                /* number of variables */
                  directions::X3 ,  /* Vary indices in X1 direction */
                  f64               /* Double precision */
                 ) - indices.numGhostX3;

    grid prim(params::N1, params::N2, params::N3,
              params::numGhost, params::dim, vars::dof,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
             );

    grid cons(params::N1, params::N2, params::N3,
              params::numGhost, params::dim, vars::dof,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
             );

    setXCoords(indices, locations::CENTER, XCoords);

    geometry geomCenter(XCoords);

    geometry geomTile(geomCenter,
                      0, params::tileSizeX1,
                      0, params::tileSizeX2,
                      0, params::tileSizeX3
                     );

    array primTile[vars::dof], consTile[vars::dof];
    af::seq tileDomainX1(0, params::tileSizeX1-1);
    af::seq tileDomainX2(0, params::tileSizeX2-1);
    af::seq tileDomainX3(0, params::tileSizeX3-1);

    /* Allocate space for tile */
    for (int var=0; var<vars::dof; var++)
    {
      primTile[var] = af::constant(0., 
                                   params::tileSizeX1,
                                   params::tileSizeX2,
                                   params::tileSizeX3, f64
                                  );
      consTile[var] = af::constant(0.,
                                   params::tileSizeX1,
                                   params::tileSizeX2,
                                   params::tileSizeX3, f64
                                  );
    }

    fluidElement elemTile(primTile, geomTile);

    int numEvals = 100;
    tiledComputation(prim, geomCenter, cons, 
                     elemTile,
                     geomTile, primTile, consTile
                    );
    af::sync();

    double numReads = 138;
    double numWrites = 38;
    double memoryBandwidth = 0.;
    double timeElapsed = 0.;

    af::timer::start();
    for (int n=0; n<numEvals; n++)
    {
      tiledComputation(prim, geomCenter, cons, 
                       elemTile, geomTile, primTile,
                       consTile
                      );
    }
    timeElapsed = af::timer::stop();
    printf("Time taken for %d tiled computation = %g secs\n", 
           numEvals, timeElapsed
          );
    memoryBandwidth =   (double)(params::N1) 
                      * (double)(params::N2)
                      * (double)(params::N3)
                      * 8. * (numReads + numWrites)/timeElapsed/1e9;

    printf("Memory bandwidth for tiled computation = %g GB/sec\n",
           memoryBandwidth
          );

    fluidElement elem(prim.vars, geomCenter);
    af::sync();

    af::timer::start();
    for (int n=0; n<numEvals; n++)
    {
      elem.set(prim.vars, geomCenter);
    }
    af::sync();
    timeElapsed = af::timer::stop();
    printf("Time taken for %d untiled computation = %g secs\n",
           numEvals, timeElapsed
          );
    memoryBandwidth =   (double)(params::N1) 
                      * (double)(params::N2)
                      * (double)(params::N3)
                      * 8. * (numReads + numWrites)/timeElapsed/1e9;

    printf("Memory bandwidth for untiled computation = %g GB/sec\n",
           memoryBandwidth
          );
    
//    timeStepper ts;
//
//    params::Time = 0.;
//    while(params::Time < params::finalTime)
//    {
//    	ts.timeStep();
//    	params::Time+=params::dt;
//     }
  }

  PetscFinalize();  
  return(0);
}

void tiledComputation(const grid &prim, 
                      const geometry &geom, 
                      grid &cons,
                      fluidElement &elemTile,
                      geometry &geomTile,
                      array primTile[vars::dof],
                      array consTile[vars::dof]
                     )
{
  for (int kTile=0; kTile<params::N3/params::tileSizeX3; kTile++)
  {
    //printf("kTile = %d of %d\n", kTile, params::N3/tileSizeX3);
    for (int jTile=0; jTile<params::N2/params::tileSizeX2; jTile++)
    {
      for (int iTile=0; iTile<params::N1/params::tileSizeX1; iTile++)
      {
        int iStart = iTile     * params::tileSizeX1;
        int iEnd   = (iTile+1) * params::tileSizeX1;

        int jStart = jTile     * params::tileSizeX2;
        int jEnd   = (jTile+1) * params::tileSizeX2;

        int kStart = kTile     * params::tileSizeX3;
        int kEnd   = (kTile+1) * params::tileSizeX3;

        af::seq tileDomainX1(iStart, iEnd-1);
        af::seq tileDomainX2(jStart, jEnd-1);
        af::seq tileDomainX3(kStart, kEnd-1);

        geomTile.copyFrom(geom, 
                          iStart, iEnd,
                          jStart, jEnd,
                          kStart, kEnd
                         );

        /* Copy data to tile */
        for (int var=0; var<vars::dof; var++)
        {
          primTile[var] = 
           prim.vars[var](tileDomainX1, tileDomainX2, tileDomainX3).copy();
        }

        elemTile.set(primTile, geomTile);
        elemTile.computeFluxes(geomTile, 0, consTile);

        /* Copy computed data from tile */
        for (int var=0; var<vars::dof; var++)
        {
          cons.vars[var](tileDomainX1, tileDomainX2, tileDomainX3) = 
            consTile[var].copy();
        }
      }
    }
  }
  af::sync();
}
