#include "grim.hpp"
#include "params.cpp"

void initialConditions(const array xCoords[3],
                       array prim[vars::dof]
                      );

/* Returns memory bandwidth in GB/sec */
double memoryBandwidth(const double numReads,
                       const double numWrites,
                       const double numEvals,
                       const double timeElapsed
                      )
{
  switch (params::dim)
  {
    case 1:
    return   (double)(params::N1 + 2*params::numGhost) 
           * 8. * (numReads + numWrites) * numEvals /timeElapsed/1e9;

    case 2:
    return   (double)(params::N1 + 2*params::numGhost) 
           * (double)(params::N2 + 2*params::numGhost)
           * 8. * (numReads + numWrites) * numEvals /timeElapsed/1e9;

    case 3:
    return   (double)(params::N1 + 2*params::numGhost) 
           * (double)(params::N2 + 2*params::numGhost)
           * (double)(params::N3 + 2*params::numGhost)
           * 8. * (numReads + numWrites) * numEvals /timeElapsed/1e9;

  }
}

void riemannSolver(fluidElement &elemFace,
                   const grid &primLeft, const grid &primRight,
                   const grid &fluxLeft, const grid &fluxRight,
                   const grid &consLeft, const grid &consRight,
                   const geometry &geomLeft, const geometry &geomRight,
                   const int dir,
                   grid &flux,
                   int &numReads, int &numWrites
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

    grid primOld(params::N1, params::N2, params::N3,
                 params::numGhost, params::dim, vars::dof,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                );

    grid primLeft(params::N1, params::N2, params::N3,
                  params::numGhost, params::dim, vars::dof,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                 );

    grid primRight(params::N1, params::N2, params::N3,
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

    grid consOld(params::N1, params::N2, params::N3,
                 params::numGhost, params::dim, vars::dof,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                );

    grid sourcesExplicit(params::N1, params::N2, params::N3,
                         params::numGhost, params::dim, vars::dof,
                         DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                         DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                         DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                        );

    grid sourcesImplicitOld(params::N1, params::N2, params::N3,
                            params::numGhost, params::dim, vars::dof,
                            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                           );

    grid sourcesImplicit(params::N1, params::N2, params::N3,
                         params::numGhost, params::dim, vars::dof,
                         DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                         DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                         DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                        );

    grid sourcesTimeDer(params::N1, params::N2, params::N3,
                        params::numGhost, params::dim, vars::dof,
                        DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                        DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                        DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                       );

    grid fluxLeft(params::N1, params::N2, params::N3,
                  params::numGhost, params::dim, vars::dof,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                 );

    grid fluxRight(params::N1, params::N2, params::N3,
                   params::numGhost, params::dim, vars::dof,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                  );

    grid consLeft(params::N1, params::N2, params::N3,
                  params::numGhost, params::dim, vars::dof,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                  DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                 );

    grid consRight(params::N1, params::N2, params::N3,
                   params::numGhost, params::dim, vars::dof,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                  );

    grid fluxX1(params::N1, params::N2, params::N3,
                params::numGhost, params::dim, vars::dof,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
               );

    grid fluxX2(params::N1, params::N2, params::N3,
                params::numGhost, params::dim, vars::dof,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
               );

    grid fluxX3(params::N1, params::N2, params::N3,
                params::numGhost, params::dim, vars::dof,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
               );

    setXCoords(indices, locations::CENTER, XCoords);
    geometry geomCenter(XCoords);
    geomCenter.computeConnectionCoeffs();

    initialConditions(XCoords.vars, prim.vars);
    initialConditions(XCoords.vars, primOld.vars);

    setXCoords(indices, locations::LEFT, XCoords);
    geometry geomLeft(XCoords);

    setXCoords(indices, locations::RIGHT, XCoords);
    geometry geomRight(XCoords);

    setXCoords(indices, locations::TOP, XCoords);
    geometry geomTop(XCoords);

    setXCoords(indices, locations::BOTTOM, XCoords);
    geometry geomBottom(XCoords);

    int numEvals = 10;
    double timeElapsed = 0.;
    int numReads, numWrites;
    int numReadsElemSet, numWritesElemSet;
    int numReadsComputeFluxes, numWritesComputeFluxes;

    fluidElement elem(prim.vars, geomCenter, 
                      numReadsElemSet, numWritesElemSet
                     );
    fluidElement elemOld(primOld.vars, geomCenter, 
                         numReadsElemSet, numWritesElemSet
                        );
    fluidElement elemFace(prim.vars, geomLeft, 
                          numReadsElemSet, numWritesElemSet
                         );
    elem.computeFluxes(geomCenter, 0, cons.vars,
                       numReadsComputeFluxes, numWritesComputeFluxes
                      );
    double dX[3];
    dX[1] = prim.dX1;
    dX[2] = prim.dX2;
    dX[3] = prim.dX3;
    for (int n=0; n<numEvals; n++)
    {
      elemOld.computeEMHDGradients(geomCenter, dX, numReads, numWrites);
      elemOld.computeExplicitSources(geomCenter, 
                                     sourcesExplicit.vars,
                                     numReads, numWrites
                                    );
      elem.computeImplicitSources(geomCenter, 
                                  sourcesImplicit.vars,
                                  numReads, 
                                  numWrites
                                 );
      elem.computeTimeDerivSources(geomCenter,
                                   elemOld, elem, params::dt,
                                   sourcesImplicit.vars,
                                   numReads, 
                                   numWrites
                                  );
      reconstruction::reconstruct(prim, directions::X1, primLeft, primRight,
                                  numReads, numWrites
                                 );
      riemannSolver(elemFace,
                    primLeft, primRight,
                    fluxLeft, fluxRight,
                    consLeft, consRight,
                    geomLeft, geomRight,
                    directions::X1,
                    fluxX1,
                    numReads, numWrites
                  );
    }
    af::sync();

    af::timer::start();
    for (int n=0; n<100*numEvals; n++)
    {
      elem.set(prim.vars, geomCenter, 
               numReadsElemSet, numWritesElemSet
              );
      elem.computeFluxes(geomCenter, 0, cons.vars,
                         numReadsComputeFluxes, numWritesComputeFluxes
                        );
    }
    af::sync();
    timeElapsed = af::timer::stop();
    numReads = numReadsElemSet + numReadsComputeFluxes;
    numWrites = numWritesElemSet + numWritesComputeFluxes;
    printf("\nConserved vars computation:\n");
    printf("Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
           100*numEvals, timeElapsed, 
           memoryBandwidth(numReads, numWrites, 100*numEvals, timeElapsed)
          );

    af::timer::start();
    for (int n=0; n<100*numEvals; n++)
    {
      elem.computeImplicitSources(geomCenter, 
                                  sourcesImplicit.vars,
                                  numReads, 
                                  numWrites
                                 );
      elemOld.computeImplicitSources(geomCenter, 
                                     sourcesImplicitOld.vars,
                                     numReads, 
                                     numWrites
                                    );
    }
    af::sync();
    timeElapsed = af::timer::stop();
    printf("\nImplicit sources computation:\n");
    printf("Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
           100*numEvals, timeElapsed, 
           memoryBandwidth(2*numReads, 2*numWrites, 100*numEvals, timeElapsed)
          );

    af::timer::start();
    for (int n=0; n<100*numEvals; n++)
    {
      elem.computeTimeDerivSources(geomCenter,
                                   elemOld, elem, params::dt,
                                   sourcesImplicit.vars,
                                   numReads, 
                                   numWrites
                                  );
    }
    af::sync();
    timeElapsed = af::timer::stop();
    printf("\nTime derivative sources computation:\n");
    printf("Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
           100*numEvals, timeElapsed, 
           memoryBandwidth(numReads, numWrites, 100*numEvals, timeElapsed)
          );

    af::timer::start();
    for (int n=0; n<numEvals; n++)
    {
      elemOld.computeExplicitSources(geomCenter, 
                                     sourcesExplicit.vars,
                                     numReads, numWrites
                                    );
    }
    af::sync();
    timeElapsed = af::timer::stop();
    printf("\nExplicit sources computation:\n");
    printf("Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
           numEvals, timeElapsed, 
           memoryBandwidth(numReads, numWrites, numEvals, timeElapsed)
          );

    for (int dir=directions::X1; dir<=directions::X3; dir++)
    {
      af::timer::start();
      for (int n=0; n<numEvals; n++)
      {
        reconstruction::reconstruct(prim, directions::X1, primLeft, primRight,
                                    numReads, numWrites
                                   );
      }
      af::sync();
      timeElapsed = af::timer::stop();
      printf("\nReconstruction dir %d\n", dir);
      printf("Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
             numEvals, timeElapsed, 
             memoryBandwidth(numReads, numWrites, numEvals, timeElapsed)
            );

      switch (dir)
      {
        case directions::X1:
          af::timer::start();
          for (int n=0; n<numEvals; n++)
          {
            riemannSolver(elemFace,
                          primLeft, primRight,
                          fluxLeft, fluxRight,
                          consLeft, consRight,
                          geomLeft, geomRight,
                          dir,
                          fluxX1,
                          numReads, numWrites
                         );
          }
          af::sync();
          timeElapsed = af::timer::stop();
          break;

        case directions::X2:
          af::timer::start();
          for (int n=0; n<numEvals; n++)
          {
            riemannSolver(elemFace,
                          primLeft, primRight,
                          fluxLeft, fluxRight,
                          consLeft, consRight,
                          geomLeft, geomRight,
                          dir,
                          fluxX2,
                          numReads, numWrites
                         );
          }
          af::sync();
          timeElapsed = af::timer::stop();
          break;

        case directions::X3:
          af::timer::start();
          for (int n=0; n<numEvals; n++)
          {
            riemannSolver(elemFace,
                          primLeft, primRight,
                          fluxLeft, fluxRight,
                          consLeft, consRight,
                          geomLeft, geomRight,
                          dir,
                          fluxX3,
                          numReads, numWrites
                         );
          }
          af::sync();
          timeElapsed = af::timer::stop();
          break;
      }

      printf("\nRiemann problem in dir %d:\n", dir);
      printf("Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
             numEvals, timeElapsed, 
             memoryBandwidth(numReads, numWrites, numEvals, timeElapsed)
            );

    }

    af::timer::start();
    for (int n=0; n<numEvals; n++)
    {
      elemOld.computeEMHDGradients(geomCenter, dX, numReads, numWrites);
    }
    af::sync();
    timeElapsed = af::timer::stop();
    printf("\nEMHD gradient computation :\n");
    printf("Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
           numEvals, timeElapsed, 
           memoryBandwidth(numReads, numWrites, numEvals, timeElapsed)
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

void riemannSolver(fluidElement &elemFace,
                   const grid &primLeft,
                   const grid &primRight,
                   const grid &fluxLeft,
                   const grid &fluxRight,
                   const grid &consLeft,
                   const grid &consRight,
                   const geometry &geomLeft,
                   const geometry &geomRight,
                   const int dir,
                   grid &flux,
                   int &numReads,
                   int &numWrites
                  )
{
  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;

  elemFace.set(primLeft.vars, geomLeft, 
               numReadsElemSet, numWritesElemSet
              );
  elemFace.computeFluxes(geomLeft, 1, fluxLeft.vars,
                         numReadsComputeFluxes, 
                         numWritesComputeFluxes
                        );
  elemFace.computeFluxes(geomLeft, 0, consLeft.vars,
                         numReadsComputeFluxes, 
                         numWritesComputeFluxes
                        );

  elemFace.set(primRight.vars, geomRight, 
               numReadsElemSet, numWritesElemSet
              );
  elemFace.computeFluxes(geomRight, 1, fluxRight.vars,
                         numReadsComputeFluxes, 
                         numWritesComputeFluxes
                        );
  elemFace.computeFluxes(geomRight, 0, consRight.vars,
                         numReadsComputeFluxes, 
                         numWritesComputeFluxes
                        );


  numReads = 2*(numReadsElemSet + 2*numReadsComputeFluxes);
  numWrites = 2*(numWritesElemSet + 2*numWritesComputeFluxes);

  int shiftX1, shiftX2, shiftX3;
  switch (dir)
  {
    case directions::X1:
      shiftX1  = 1;
      shiftX2  = 0;
      shiftX3  = 0;
      break;

    case directions::X2:
      shiftX1  = 0;
      shiftX2  = 1;
      shiftX3  = 0;
      break;

    case directions::X3:
      shiftX1  = 0;
      shiftX2  = 0;
      shiftX3  = 1;
      break;
  }

  for (int var=0; var<vars::dof; var++)
  {
    flux.vars[var] = 
      af::shift(fluxLeft.vars[var], shiftX1, shiftX2, shiftX3)
    - fluxRight.vars[var]
    + consRight.vars[var] 
    - af::shift(consLeft.vars[var], shiftX1, shiftX2, shiftX3);

    flux.vars[var].eval();
  }
  /* Reads:
   * -----
   *  fluxLeft[var], fluxRight[var], consLeft[var], consRight[var] : 4*vars::dof
   *
   * Writes:
   * ------
   * flux[var] : vars::dof */
  numReads  += 4*vars::dof;
  numWrites +=  vars::dof;
}

void initialConditions(const array xCoords[3],
                       array prim[vars::dof]
                      )
{
  double Aw = 1.e-5;
  double k1 = 2.*M_PI;
  double k2 = 4.*M_PI;
  double Gamma = - 0.5533585207638141;
  double Omega = - 3.6262571286888425;

  array cphi = af::cos(  k1*xCoords[directions::X1]
                  		 + k2*xCoords[directions::X2]
                      );

  array sphi = af::sin(  k1*xCoords[directions::X1]
                  		 + k2*xCoords[directions::X2]
                      );

  /* Initial conditions */

  //Full EMHD mode (from grim2D)
  prim[vars::RHO] = 1.;
  prim[vars::U]   = 2.;
  prim[vars::U1]  = 0.;
  prim[vars::U2]  = 0.; 
  prim[vars::U3]  = 0.; 
  prim[vars::B1]  = 0.1; 
  prim[vars::B2]  = 0.3;
  prim[vars::B3]  = 0.;
  prim[vars::Q]  = 0.;
  prim[vars::DP]  = 0.;

  prim[vars::RHO] += Aw*cphi*(-0.518522524082246)
                    +Aw*sphi*0.1792647678001878;

  prim[vars::U]   += Aw*cphi*0.5516170736393813;

  prim[vars::U1]  += Aw*cphi*0.008463122479547856
                    +Aw*sphi*(-0.011862022608466367);

  prim[vars::U2]  += Aw*cphi*(-0.16175466371870734)
                    +Aw*sphi*(0.034828080823603294);

  prim[vars::B1]  += Aw*cphi*(-0.05973794979640743)
                    +Aw*sphi*0.03351707506150924;

  prim[vars::B2]  += Aw*cphi*0.02986897489820372
                    -Aw*sphi*0.016758537530754618;

  prim[vars::Q]   += Aw*cphi*0.5233486841539436
                    -Aw*sphi*0.04767672501939603;

  prim[vars::DP]  += Aw*cphi*0.2909106062057657
                    -Aw*sphi*0.02159452055336572;

  for (int var=0; var<vars::dof; var++)
  {
    prim[var].eval();
  }

  af::sync();
}
