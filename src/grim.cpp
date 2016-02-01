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

void bandwidthTest()
{
//  array a = af::randu(params::N1, params::N2, params::N3);
//  array b = af::randu(params::N1, params::N2, params::N3);
//  array c = af::constant(0., params::N1, params::N2, params::N3);
//  a.eval(); b.eval(); c.eval();
//  af::sync();

  grid prim(params::N1, params::N2, params::N3,
            params::numGhost, params::dim, 3,
            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
           );

  int numEvals = 10000;
  af::timer::start();
  for (int n=0; n < numEvals; n++)
  {
    prim.vars[2] = prim.vars[1] + prim.vars[0];
    prim.vars[2].eval();
  }
  af::sync();
  double timeElapsed = af::timer::stop();
  PetscPrintf(PETSC_COMM_WORLD, "Num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
         numEvals, timeElapsed, 
         memoryBandwidth(2, 1, numEvals, timeElapsed)
        );
}

void computeResidual(const grid &prim,
                     fluidElement *elem,
                     fluidElement *elemOld,
                     fluidElement *elemHalfStep,
                     geometry *geom,
                     grid *sourcesExplicit,
                     grid *sourcesImplicit,
                     grid *sourcesImplicitOld,
                     grid *sourcesTimeDer,
                     grid *cons,
                     grid *consOld,
                     grid *divFluxes,
                     grid &residual, 
                     const bool ComputeExplicitTerms,
                     const int currentStep
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
    //af::setDevice(4);
    af::info();
  };

  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    grid XCoords(params::N1, params::N2, params::N3,
                 params::numGhost, params::dim, 3,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                );

    grid prim(params::N1, params::N2, params::N3,
              params::numGhost, params::dim, vars::dof,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
              DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
             );

    grid residual(params::N1, params::N2, params::N3,
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

    grid divFluxes(params::N1, params::N2, params::N3,
                   params::numGhost, params::dim, vars::dof,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                   DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
                  );

    setXCoords(locations::CENTER, XCoords);
    geometry geomCenter(XCoords);
    geomCenter.computeConnectionCoeffs();

    initialConditions(XCoords.vars, prim.vars);
    initialConditions(XCoords.vars, primOld.vars);

    setXCoords(locations::LEFT, XCoords);
    geometry geomLeft(XCoords);

    setXCoords(locations::RIGHT, XCoords);
    geometry geomRight(XCoords);

    setXCoords(locations::TOP, XCoords);
    geometry geomTop(XCoords);

    setXCoords(locations::BOTTOM, XCoords);
    geometry geomBottom(XCoords);

    riemannSolver riemann(prim, geomCenter);

    int numEvals = 10;
    double timeElapsed = 0.;
    int numReads, numWrites;
    int numReadsElemSet, numWritesElemSet;
    int numReadsComputeFluxes, numWritesComputeFluxes;

    fluidElement elem(prim.vars, geomCenter, 
                      numReadsElemSet, numWritesElemSet
                     );
    fluidElement elemHalfStep(prim.vars, geomCenter, 
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
    for (int n=0; n<1; n++)
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
      riemann.solve(primLeft, primRight,
                    geomLeft, geomRight,
                    directions::X1,
                    fluxX1,
                    numReads, numWrites
                   );
      for (int var=0; var<vars::dof; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxX1.dX1);
        array filterX2 = array(1, 3, 1, 1, filter1D)/(fluxX2.dX2);
        array filterX3 = array(1, 1, 3, 1, filter1D)/(fluxX3.dX3);

        array dFluxX1_dX1 = convolve(fluxX1.vars[var], filterX1);
        array dFluxX2_dX2 = convolve(fluxX2.vars[var], filterX2);
        array dFluxX3_dX3 = convolve(fluxX3.vars[var], filterX3);

        divFluxes.vars[var] = dFluxX1_dX1 + dFluxX2_dX2 + dFluxX3_dX3;
        divFluxes.vars[var].eval();
      }
      computeResidual(prim, &elem, &elemOld, &elemHalfStep,
                      &geomCenter, &sourcesExplicit, &sourcesImplicit,
                      &sourcesImplicitOld, &sourcesTimeDer,
                      &cons, &consOld, &divFluxes, residual, true, 1
                     );
    }
    af::sync();

    printf("\nKernel compilation complete\n");


    printf("\n Performing bandwidth test...\n");
    bandwidthTest();

    af::timer::start();
    for (int n=0; n<numEvals; n++)
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
           numEvals, timeElapsed, 
           memoryBandwidth(numReads, numWrites, numEvals, timeElapsed)
          );

    af::timer::start();
    for (int n=0; n<numEvals; n++)
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
           numEvals, timeElapsed, 
           memoryBandwidth(2*numReads, 2*numWrites, numEvals, timeElapsed)
          );

    af::timer::start();
    for (int n=0; n<numEvals; n++)
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
           numEvals, timeElapsed, 
           memoryBandwidth(numReads, numWrites, numEvals, timeElapsed)
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
            riemann.solve(primLeft, primRight,
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
            riemann.solve(primLeft, primRight,
                          geomBottom, geomTop,
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
            riemann.solve(primLeft, primRight,
                          geomCenter, geomCenter,
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

    af::timer::start();
    for (int n=0; n<numEvals; n++)
    {
      for (int var=0; var<vars::dof; var++)
      {
        double filter1D[] = {1, -1, 0}; /* Forward difference */
    
        array filterX1 = array(3, 1, 1, 1, filter1D)/(fluxX1.dX1);
        array filterX2 = array(1, 3, 1, 1, filter1D)/(fluxX2.dX2);
        array filterX3 = array(1, 1, 3, 1, filter1D)/(fluxX3.dX3);

        array dFluxX1_dX1 = convolve(fluxX1.vars[var], filterX1);
        array dFluxX2_dX2 = convolve(fluxX2.vars[var], filterX2);
        array dFluxX3_dX3 = convolve(fluxX3.vars[var], filterX3);

        divFluxes.vars[var] = dFluxX1_dX1 + dFluxX2_dX2 + dFluxX3_dX3;
        divFluxes.vars[var].eval();
      }
    }
    af::sync();
    timeElapsed = af::timer::stop();
    printf("\ndivFluxes computation :\n");
    printf("Num evals = %d, time taken = %g secs\n", numEvals, timeElapsed);

    af::timer::start();
    for (int n=0; n<numEvals; n++)
    {
      computeResidual(prim, &elem, &elemOld, &elemHalfStep,
                      &geomCenter, &sourcesExplicit, &sourcesImplicit,
                      &sourcesImplicitOld, &sourcesTimeDer,
                      &cons, &consOld, &divFluxes, residual, true, 1
                     );
    }
    af::sync();
    timeElapsed = af::timer::stop();
    printf("\nResidual computation :\n");
    printf("Num evals = %d, time taken = %g secs\n", numEvals, timeElapsed);

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

void computeResidual(const grid &prim,
                     fluidElement *elem,
                     fluidElement *elemOld,
                     fluidElement *elemHalfStep,
                     geometry *geom,
                     grid *sourcesExplicit,
                     grid *sourcesImplicit,
                     grid *sourcesImplicitOld,
                     grid *sourcesTimeDer,
                     grid *cons,
                     grid *consOld,
                     grid *divFluxes,
                     grid &residual, 
                     const bool ComputeExplicitTerms,
                     const int currentStep
                    )
{
  int numReads, numWrites;
  elem->set(prim.vars, *geom, numReads, numWrites);
  elem->computeFluxes(*geom, 0, cons->vars, numReads, numWrites);

  if (currentStep == 1)
  {
    int useImplicitSources = 0;
    if (ComputeExplicitTerms)
    {
      elemOld->computeExplicitSources(*geom, sourcesExplicit->vars,
                                      numReads, numWrites                            
                                     );
    }
    else
    {
      for (int var=0; var<vars::dof; var++)
      {
	      sourcesExplicit->vars[var]=0.;
      }
    }

    elemOld->computeImplicitSources(*geom, sourcesImplicitOld->vars,
                                    numReads, numWrites
                                   );
    elem->computeImplicitSources(*geom, sourcesImplicit->vars,
                                 numReads, numWrites
                                );
    elemOld->computeTimeDerivSources(*geom,
				                             *elemOld, *elem,
		                         		     params::dt/2,
                        				     sourcesTimeDer->vars,
                                     numReads, numWrites
                                    );

    for (int var=0; var<vars::dof; var++)
    {
      residual.vars[var] = 
        (cons->vars[var] - consOld->vars[var])/(params::dt/2.)
  	  + divFluxes->vars[var]
      + sourcesExplicit->vars[var]
  	  + 0.5*(sourcesImplicitOld->vars[var] + sourcesImplicit->vars[var])
	    + sourcesTimeDer->vars[var];
    }

    //Normalization of the residual
    if (params::conduction)
    {
	    if(params::highOrderTermsConduction)
      {
	      residual.vars[vars::Q] *=
          elemOld->temperature 
        * af::sqrt(elemOld->rho*elemOld->chi_emhd*elemOld->tau);
      }
	    else
      {
	      residual.vars[vars::Q] *= elemOld->tau;
      }
    }

    if (params::viscosity)
    {
	    if(params::highOrderTermsViscosity)
      {
	      residual.vars[vars::DP] *=
          af::sqrt(   elemOld->rho*elemOld->nu_emhd
                    * elemOld->temperature*elemOld->tau
                  );
      }
	    else
      {
	      residual.vars[vars::DP] *= elemOld->tau;
      }
    }

  } /* End of timeStepperSwitches::HALF_STEP */

  else if (currentStep == 2)
  {
    int useImplicitSources = 0;
    if(ComputeExplicitTerms)
    {
      elemHalfStep->computeExplicitSources(*geom, sourcesExplicit->vars,
                                           numReads, numWrites
                                          );
    }
    else
    {
      for (int var=0; var<vars::dof; var++)
      {
	      sourcesExplicit->vars[var]=0.;
      }
    }

    elemOld->computeImplicitSources(*geom, sourcesImplicitOld->vars,
                                    numReads, numWrites
                                   );
    elem->computeImplicitSources(*geom, sourcesImplicit->vars,
                                 numReads, numWrites
                                );
    elemHalfStep->computeTimeDerivSources(*geom,
	                                			  *elemOld, *elem,
					                                params::dt,
					                                sourcesTimeDer->vars,
                                          numReads, numWrites
                                         );

    for (int var=0; var<vars::dof; var++)
    {
      residual.vars[var] = 
        (cons->vars[var] - consOld->vars[var])/params::dt
    	+ divFluxes->vars[var]
	    + sourcesExplicit->vars[var]
	    + 0.5*(sourcesImplicitOld->vars[var] + sourcesImplicit->vars[var])
	    + sourcesTimeDer->vars[var];
    }

    //Normalization of the residual
    if (params::conduction)
    {
	    if(params::highOrderTermsConduction)
      {
	      residual.vars[vars::Q] *= 
          elemHalfStep->temperature
        * af::sqrt(elemHalfStep->rho*elemHalfStep->chi_emhd*elemHalfStep->tau);
      }
    	else
      {
	      residual.vars[vars::Q] *= elemHalfStep->tau;
      }
    }

    if (params::viscosity)
    {
      if (params::highOrderTermsViscosity)
      {
	      residual.vars[vars::DP] *= 
          af::sqrt(   elemHalfStep->rho*elemHalfStep->nu_emhd
                    * elemHalfStep->temperature*elemHalfStep->tau
                  );
      }
	    else
      {
	      residual.vars[vars::DP] *= elemHalfStep->tau;
      }
    }

  } /* End of timeStepperSwitches::FULL_STEP */

  for (int var=0; var<vars::dof; var++)
  {
    residual.vars[var].eval();
  }
}
