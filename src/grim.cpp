#include "grim.hpp"
#include "params.hpp"

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);

  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  af::setDevice(world_rank%params::numDevices);

  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    //grid testvar(params::N1, params::N2, params::N3, params::dim,
    //             3, params::numGhost, 0, 0, 0);
    
    //testvar.vars[0] = 1.;
    //testvar.vars[1] = 2.;
    //testvar.vars[2] = testvar.vars[0] + testvar.vars[1];
    //testvar.vars[2] += 3;

    //double testsum = af::sum<double>(testvar.vars[2]);    
    //af_print(af::sum(af::sum(testvar.vars[2])));
    //PetscPrintf(PETSC_COMM_WORLD, "testsum = %e\n", testsum);
    
    
    timeStepper ts(params::N1, params::N2, params::N3,
                   params::dim, vars::dof, params::numGhost,
                   params::Time, params::InitialDt,
                   params::boundaryLeft,  params::boundaryRight,
                   params::boundaryTop,   params::boundaryBottom,
                   params::boundaryFront, params::boundaryBack,
                   params::metric, params::blackHoleSpin, params::hSlope,
                   params::X1Start, params::X1End,
                   params::X2Start, params::X2End,
                   params::X3Start, params::X3End
                  );

    
    int numReads, numWrites;
    fluidElement elem(*ts.primOld, *ts.geomCenter, numReads, numWrites);
    
    array tstCon[NDIM];
    elem.constructTetrads(*ts.geomCenter);
    elem.coordConToTetradCon(elem.uCon, tstCon);
    
    af_print(elem.uCon[0]);
    af_print(tstCon[0]);
    af_print(tstCon[1]);
    af_print(tstCon[2]);
    af_print(tstCon[3]);
    
    //elem.set()
    
    
    /*PetscPrintf(PETSC_COMM_WORLD, "  Generating compute kernels...\n\n");
    int numReads, numWrites;
    ts.timeStep(numReads, numWrites);

    af::sync();

    PetscPrintf(PETSC_COMM_WORLD, "\n  Kernel compilation complete\n");

    int n=0;
    while(ts.time<params::finalTime)
    {
      n++;
      PetscPrintf(PETSC_COMM_WORLD, "\n|----Time step %d----|\n", n);
      ts.timeStep(numReads, numWrites);
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n===Program execution complete===\n\n");
    }*/
  }

  PetscFinalize();  
  return(0);
}
