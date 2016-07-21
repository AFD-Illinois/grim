#include "../problem.hpp"

namespace problemParams
{
  std::string rhoInputFile      = "atmosphere_soln_rho.txt";
  std::string uInputFile        = "atmosphere_soln_u.txt"  ;
  std::string qInputFile        = "atmosphere_soln_phi.txt";
  std::string rCoordsInputFile  = "atmosphere_soln_rCoords.txt";
};

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  array xCoords[3],XCoords[3];
  geom.getXCoords(XCoords);
  geom.XCoordsToxCoords(XCoords,xCoords);

  array r = xCoords[directions::X1];

  chi_emhd = af::pow(r, 0.5);
  tau      = 1.2 * chi_emhd/ (soundSpeed * soundSpeed);
}

void timeStepper::initialConditions(int &numReads, int &numWrites)
{
  PetscPrintf(PETSC_COMM_WORLD, "Setting up atmosphere...");

  int numGhost = params::numGhost;
  double rho1D[N1+2*numGhost];
  double u1D[N1+2*numGhost];
  double q1D[N1+2*numGhost];
  double rCoords1D[N1+2*numGhost];

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); /* get current proc id */

  if (rank==0)
  {
    FILE *rhoFile;
    FILE *uFile;
    FILE *qFile;
    FILE *rCoordsFile;
  
    char *rhoLine     = NULL;
    char *uLine       = NULL;
    char *qLine     = NULL; 
    char *rCoordsLine = NULL;

    size_t rhoLen     = 0; ssize_t rhoRead;
    size_t uLen       = 0; ssize_t uRead;
    size_t qLen       = 0; ssize_t qRead;
    size_t rCoordsLen = 0; ssize_t rCoordsRead;

    rhoFile     = fopen(problemParams::rhoInputFile.c_str(), "r");
    uFile       = fopen(problemParams::uInputFile.c_str()  , "r");
    qFile       = fopen(problemParams::qInputFile.c_str(), "r");
    rCoordsFile = fopen(problemParams::rCoordsInputFile.c_str(), "r");

    if (   rhoFile      == NULL 
        || uFile        == NULL
        || qFile        == NULL
        || rCoordsFile  == NULL
       )
    {
      PetscPrintf(PETSC_COMM_WORLD, "Input data files not found!\n");
      exit(1);
    }

    for (int i=-numGhost; i<N1+numGhost; i++)
    {
      rhoRead     = getline(&rhoLine    , &rhoLen     , rhoFile);
      uRead       = getline(&uLine      , &uLen       , uFile);
      qRead       = getline(&qLine      , &qLen       , qFile);
      rCoordsRead = getline(&rCoordsLine, &rCoordsLen , rCoordsFile);

      if (   rhoRead      == -1
          || uRead        == -1
          || qRead        == -1
          || rCoordsRead  == -1
         )
      {
        PetscPrintf(PETSC_COMM_WORLD,
        "Found the input data files but error in reading them! Check number of grid zones in python script\n");
        exit(1);
      }

      rho1D[i+numGhost]     = atof(rhoLine);
      u1D[i+numGhost]       = atof(uLine);
      q1D[i+numGhost]       = atof(qLine);
      rCoords1D[i+numGhost] = atof(rCoordsLine);
    }

    free(rhoLine);
    free(uLine);
    free(qLine);
    free(rCoordsLine);

    fclose(rhoFile);
    fclose(uFile);
    fclose(qFile);
    fclose(rCoordsFile);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "file read complete\n");

  /* Broadcast the data from rank 0 proc to all other procs */
  MPI_Bcast(&rho1D[0]     , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&u1D[0]       , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&q1D[0]       , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&rCoords1D[0] , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  /* Broadcast complete */
  PetscPrintf(PETSC_COMM_WORLD, "Broadcast complete\n");

  array &rho = primOld->vars[vars::RHO];
  array &u   = primOld->vars[vars::U];
  array &u1  = primOld->vars[vars::U1];
  array &u2  = primOld->vars[vars::U2];
  array &u3  = primOld->vars[vars::U3];
  array &B1  = primOld->vars[vars::B1];
  array &B2  = primOld->vars[vars::B2];
  array &B3  = primOld->vars[vars::B3];

  array xCoords[3], XCoords[3];
  geomCenter->getXCoords(XCoords);
  geomCenter->XCoordsToxCoords(XCoords,xCoords);
  array r = xCoords[directions::X1];

  for (int i=0; i<primOld->N1Total; i++)
  {
    int iGlobal = i + primOld->iLocalStart;

    if (fabs(r(i, 0, 0).scalar<double>() - rCoords1D[iGlobal]) > 1e-10)
    {
      PetscPrintf(PETSC_COMM_WORLD,
                 "Mismatch in rCoords! Check r coords in python script. r = %.18f, rCoords = %.18f, i = %d\n",
                 r(i, 0, 0).scalar<double>(), rCoords1D[iGlobal], i
                 );
      exit(1);
    }
    rho(i, span, span) = rho1D[iGlobal];
    u(i, span, span)   = u1D[iGlobal];

    if (params::conduction)
    {
      array &q = primOld->vars[vars::Q];
      q(i, span, span) = q1D[iGlobal];
    }
  }

  array uCon0 = 1./af::sqrt(-geomCenter->gCov[0][0]);
  array uCon1 = 0.*u1;
  array uCon2 = 0.*u2;
  array uCon3 = 0.*u3;

  /* Formula to output vUpr in Modified KerrSchild from uUpr in Boyer-Lindquist 
   * from Ben Ryan */
  array a = geomCenter->gCov[1][1];
  array b = geomCenter->gCon[0][1];
  array c = geomCenter->gCon[0][0];
  array vUprMKS = 
    (  c*uCon1/r 
     - sqrt(-a*b*b*b*b
            -a*b*b*c*uCon1*uCon1/(r*r) - b*b*c
           )
    )/(a*b*b + c);

  /* Set values in both bulk and boundaries. Later set values in bulk where we
   * need to add random fluctuations */
  u1 = vUprMKS;

  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);
  int numGhostX1 = primOld->numGhostX1;
  int numGhostX2 = primOld->numGhostX2;
  int numGhostX3 = primOld->numGhostX3;
  int N1Local    = primOld->N1Local;
  int N2Local    = primOld->N2Local;
  int N3Local    = primOld->N3Local;

  for (int k=numGhostX3; k < N3Local+numGhostX3; k++)
  {
    for (int j=numGhostX2; j < N2Local+numGhostX2; j++)
    {
      for (int i=numGhostX1; i < N1Local+numGhostX1; i++)
      {
        double randNum = 0.5;
        PetscRandomGetValue(randNumGen, &randNum);

        u1(i,j,k) = 
          vUprMKS(i,j,k)
        * (1. + params::InitialPerturbationAmplitude*(randNum-0.5));
      }
    }
  }
  PetscRandomDestroy(&randNumGen);

  //B2 = 1e-5/geomCenter->g;
  B3 = 1e-5/geomCenter->g;

  primIC->vars[vars::RHO] = rho;
  primIC->vars[vars::U]   = u;
  primIC->vars[vars::U1]  = u1;
  
  PetscPrintf(PETSC_COMM_WORLD, "done\n");
  /* Done with setting the initial conditions */

  af::seq domainX1 = *primOld->domainX1;
  af::seq domainX2 = *primOld->domainX2;
  af::seq domainX3 = *primOld->domainX3;
  
  computeDivB(*primOld, numReads, numWrites);
  double divBNorm = 
    af::norm(af::flat(divB->vars[0](domainX1, domainX2, domainX3)), AF_NORM_VECTOR_1);
  printf("divB = %g\n", divBNorm);
}



void timeStepper::setProblemSpecificBCs(
                                        int &numReads, 
                                        int &numWrites
                                       )
{
  int numGhost = params::numGhost;
  // 1) Choose which primitive variables are corrected.
  grid* primBC;
  fluidElement* elemBC;
  if(currentStep == timeStepperSwitches::HALF_STEP)
    {
      primBC = primOld;
    }
  else
    {
      primBC = primHalfStep;
    }
  af::seq leftBoundary(0, numGhost);
  af::seq rightBoundary(primOld->N1Local + numGhost, 
                        primOld->N1Local + 2*numGhost-1
                       );

  af::seq bottomBoundary(0, numGhost);
  af::seq topBoundary(primOld->N2Local + numGhost,
                      primOld->N2Local + 2*numGhost-1
                     );

  if (primOld->iLocalStart == 0)
  {
    primBC->vars[vars::RHO](leftBoundary, span, span) = 
      primIC->vars[vars::RHO](leftBoundary, span, span);

    primBC->vars[vars::U](leftBoundary,   span, span) = 
      primIC->vars[vars::U](leftBoundary, span, span);

    primBC->vars[vars::U1](leftBoundary,  span, span) = 
      primIC->vars[vars::U1](leftBoundary, span, span);

    primBC->vars[vars::U2](leftBoundary,  span, span) = 0.;
    primBC->vars[vars::U3](leftBoundary,  span, span) = 0.;

    primBC->vars[vars::B1](leftBoundary, span, span)  = 
      0e-5/geomCenter->g(leftBoundary, span, span) ;
    primBC->vars[vars::B2](leftBoundary, span, span)  = 
      0e-5/geomCenter->g(leftBoundary, span, span) ;
    primBC->vars[vars::B3](leftBoundary, span, span)  = 
      0e-5/geomCenter->g(leftBoundary, span, span) ;
    
    if (params::conduction)
    {
      primBC->vars[vars::Q](leftBoundary, span, span) = 0.;
    }
  }

  if (primOld->iLocalEnd == N1)
  {
    primBC->vars[vars::RHO](rightBoundary, span, span) =
      primIC->vars[vars::RHO](rightBoundary, span, span);

    primBC->vars[vars::U](rightBoundary,   span, span) =
      primIC->vars[vars::U](rightBoundary, span, span);

    primBC->vars[vars::U1](rightBoundary,  span, span) =
      primIC->vars[vars::U1](rightBoundary, span, span);

    primBC->vars[vars::U2](rightBoundary, span, span) = 0.;
    primBC->vars[vars::U3](rightBoundary, span, span) = 0.;

    primBC->vars[vars::B1](rightBoundary, span, span)  = 
      0e-5/geomCenter->g(rightBoundary, span, span) ;
    primBC->vars[vars::B2](rightBoundary, span, span)  = 
      0e-5/geomCenter->g(rightBoundary, span, span) ;
    primBC->vars[vars::B3](rightBoundary, span, span)  = 
      0e-5/geomCenter->g(rightBoundary, span, span) ;

    if (params::conduction)
    {
      primBC->vars[vars::Q](rightBoundary, span, span) = 0.;
    }
  }

//  if (primOld->jLocalStart == 0)
//  {
//    B1Left.vars[0](span, bottomBoundary, span) = 0.;
//    B2Bottom.vars[0](span, bottomBoundary, span) =
//      0.*1e-5/geomBottom->g(span, bottomBoundary, span); ;
//    B3Back.vars[0](span, bottomBoundary, span) =
//      0e-5/geomCenter->g(span, bottomBoundary, span);
//
//    if (params::conduction)
//    {
//      prim.vars[vars::Q](rightBoundary, span, span) = 0.;
//    }
//  }
//
//  if (primOld->jLocalEnd == N2)
//  {
//    B1Left.vars[0](span, topBoundary, span) = 0.;
//    B2Bottom.vars[0](span, topBoundary, span) =
//      0.*1e-5/geomBottom->g(span, topBoundary, span);
//    B3Back.vars[0](span, topBoundary, span) =
//      0e-5/geomCenter->g(span, topBoundary, span);
//
//    if (params::conduction)
//    {
//      prim.vars[vars::Q](rightBoundary, span, span) = 0.;
//    }
//  }

}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{
}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
   bool WriteData = (floor(time/params::WriteDataEveryDt) != floor((time-dt)/params::WriteDataEveryDt));
  if(WriteData)
  {
    const int WriteIdx = floor(time/params::WriteDataEveryDt);
    if(WriteIdx==0)
	  {
	    PetscPrintf(PETSC_COMM_WORLD, "Printing gCov\n");
	    geomCenter->gCovGrid->dump("gCov","gCovCenter.h5");
	    geomCenter->gConGrid->dump("gCon","gConCenter.h5");
	    geomCenter->gGrid->dump("sqrtDetg","sqrtDetgCenter.h5");
	    geomCenter->xCoordsGrid->dump("xCoords","xCoordsCenter.h5");

	    geomLeft->gCovGrid->dump("gCov","gCovLeft.h5");
	    geomLeft->gConGrid->dump("gCon","gConLeft.h5");
	    geomLeft->gGrid->dump("sqrtDetg","sqrtDetgLeft.h5");
	    geomLeft->xCoordsGrid->dump("xCoords","xCoordsLeft.h5");

	    geomBottom->gCovGrid->dump("gCov","gCovBottom.h5");
	    geomBottom->gConGrid->dump("gCon","gConBottom.h5");
	    geomBottom->gGrid->dump("sqrtDetg","sqrtDetgBottom.h5");
	    geomBottom->xCoordsGrid->dump("xCoords","xCoordsBottom.h5");
	  }
      
    std::string filename   = "primVarsT";
    std::string filenameVTS= "primVarsT";
    std::string B1filename = "B1Left";
    std::string B2filename = "B2Bottom";
    std::string B3filename = "B3Back";
    std::string s_idx = std::to_string(WriteIdx);
      
    for(int i=0;i<6-s_idx.size();i++)
    {
	    filename=filename+"0";
	    B1filename=B1filename+"0";
	    B2filename=B2filename+"0";
	    B3filename=B3filename+"0";
    }
      
    filename=filename+s_idx;
    filenameVTS=filename;
    filename=filename+".h5";
    primOld->dump("primitives",filename);

    filenameVTS=filenameVTS+".vts";
    std::string varNames[vars::dof];
    varNames[vars::RHO] = "rho";
    varNames[vars::U]   = "u";
    varNames[vars::U1]  = "u1";
    varNames[vars::U2]  = "u2";
    varNames[vars::U3]  = "u3";
    varNames[vars::B1]  = "B1";
    varNames[vars::B2]  = "B2";
    varNames[vars::B3]  = "B3";
    if (params::conduction)
    {
      varNames[vars::Q]   = "q";
    }
    if (params::viscosity)
    {
      varNames[vars::DP]  = "dP";
    }

    primOld->dumpVTS(*geomCenter->xCoordsGrid, varNames, filenameVTS);

  }

//  computeDivB(*primOld, numReads, numWrites);
//  double divBNorm = 
//    af::norm(af::flat(divB->vars[0](domainX1, domainX2, domainX3)), AF_NORM_VECTOR_1);
//  printf("divB = %g\n", divBNorm);

}

void timeStepper::applyProblemSpecificFluxFilter(int &numReads,int &numWrites)
{

}
