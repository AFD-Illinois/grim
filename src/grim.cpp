#include "grim.hpp"
#include "params.cpp"

void tiledComputation(const grid &prim, const geometry &geom, grid &cons);

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


    tiledComputation(prim, geomCenter, cons);
    af::sync();

    int numEvals = 10;
    af::timer::start();
    for (int n=0; n<numEvals; n++)
    {
      tiledComputation(prim, geomCenter, cons);
    }
    printf("Time taken for %d tiled computation = %g secs\n", numEvals, af::timer::stop());

//    fluidElement elem(prim.vars, geomCenter);
//    elem.set(prim.vars, geomCenter);
//    elem.computeFluxes(geomCenter, 0, cons.vars);
//    af::sync();
//
//    af::timer::start();
//    for (int n=0; n<numEvals; n++)
//    {
//      elem.set(prim.vars, geomCenter);
//      elem.computeFluxes(geomCenter, 0, cons.vars);
//      af::sync();
//    }
//    printf("Time taken for %d untiled computation = %g secs\n", numEvals, af::timer::stop());

    
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

void tiledComputation(const grid &prim, const geometry &geom, grid &cons)
{
  int tileSizeX1 = 256;
  int tileSizeX2 = 256;
  int tileSizeX3 = 50;
  geometry geomTile(geom,
                    0, tileSizeX1,
                    0, tileSizeX2,
                    0, tileSizeX3
                   );

  array primTile[vars::dof], consTile[vars::dof];
  af::seq tileDomainX1(0, tileSizeX1-1);
  af::seq tileDomainX2(0, tileSizeX2-1);
  af::seq tileDomainX3(0, tileSizeX3-1);

  /* Allocate space for tile */
  for (int var=0; var<vars::dof; var++)
  {
    primTile[var] = af::constant(0., tileSizeX1, tileSizeX2, tileSizeX3, f64);
    consTile[var] = af::constant(0., tileSizeX1, tileSizeX2, tileSizeX3, f64);
  }

  /* Copy data to tile */
  for (int var=0; var<vars::dof; var++)
  {
    primTile[var] = 
      prim.vars[var](tileDomainX1, tileDomainX2, tileDomainX3).copy();
  }

  fluidElement elemTile(primTile, geomTile);

  for (int kTile=0; kTile<params::N3/tileSizeX3; kTile++)
  {
    //printf("kTile = %d of %d\n", kTile, params::N3/tileSizeX3);
    for (int jTile=0; jTile<params::N2/tileSizeX2; jTile++)
    {
      for (int iTile=0; iTile<params::N1/tileSizeX1; iTile++)
      {
        int iStart = iTile     * tileSizeX1;
        int iEnd   = (iTile+1) * tileSizeX1;

        int jStart = jTile     * tileSizeX2;
        int jEnd   = (jTile+1) * tileSizeX2;

        int kStart = kTile     * tileSizeX3;
        int kEnd   = (kTile+1) * tileSizeX3;

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

        /* Copy data to tile */
        for (int var=0; var<vars::dof; var++)
        {
          cons.vars[var](tileDomainX1, tileDomainX2, tileDomainX3) = 
            consTile[var].copy();
        }
      }
    }
    af::sync();
  }
}
