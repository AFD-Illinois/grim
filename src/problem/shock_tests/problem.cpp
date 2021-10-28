#include "../problem.hpp"

namespace problemParams
{
  std::string rhoInputFile     = "shock_soln_rho.txt";
  std::string uInputFile       = "shock_soln_u.txt"  ;
  std::string u1InputFile      = "shock_soln_u1.txt"  ;
  std::string qInputFile       = "shock_soln_q.txt";
  std::string dPInputFile      = "shock_soln_dP.txt"  ;
  std::string xCoordsInputFile = "shock_soln_xCoords.txt";
};

void fluidElement::setFluidElementParameters()
{
  tau = 0.1 * one;
  chi_emhd = params::ConductionAlpha * soundSpeed*soundSpeed*tau;
  nu_emhd  = params::ViscosityAlpha  * soundSpeed*soundSpeed*tau;
}

void timeStepper::initialConditions(int &numReads,int &numWrites)
{
  double X1Center = (params::X1Start + params::X1End)/2.;

  double rhoLeft,       rhoRight;
  double pressureLeft,  pressureRight;
  double u1Left,        u1Right;
  double u2Left,        u2Right;
  double u3Left,        u3Right;
  double B1Left,        B1Right;
  double B2Left,        B2Right;
  double B3Left,        B3Right;

  if (params::shockTest == "fast_shock")
  {
    params::finalTime = 2.5;
    rhoLeft       = 1.;         rhoRight      = 25.48;
    pressureLeft  = 1.;         pressureRight = 367.5;
    u1Left        = 25.;        u1Right       = 1.091;
    u2Left        = 0.;         u2Right       = 0.3923;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 20.;        B1Right       = 20.;
    B2Left        = 25.02;      B2Right       = 49.;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "slow_shock")
  {
    params::finalTime = 2.;
    rhoLeft       = 1.;         rhoRight      = 3.323;
    pressureLeft  = 10.;        pressureRight = 55.36;
    u1Left        = 1.53;       u1Right       = 0.9571;
    u2Left        = 0.;         u2Right       = -0.6822;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 10.;        B1Right       = 10.;
    B2Left        = 18.28;      B2Right       = 14.49;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "switch_on_slow")
  {
    params::finalTime = 2.;
    rhoLeft       = 1.78e-3;    rhoRight      = 0.01;
    pressureLeft  = .1;         pressureRight = 1.;
    u1Left        = -.765;      u1Right       = 0.;
    u2Left        = -1.386;     u2Right       = 0.;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 1.;         B1Right       = 1.;
    B2Left        = 1.022;      B2Right       = 0.;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "switch_off_fast")
  {
    params::finalTime = 1.;
    rhoLeft       = 0.1;        rhoRight      = 0.562;
    pressureLeft  = 1.;         pressureRight = 10.;
    u1Left        = -2.;        u1Right       = -0.212;
    u2Left        = 0.;         u2Right       = -0.590;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 2.;         B1Right       = 2.;
    B2Left        = 0.;         B2Right       = 4.71;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "alfven_wave")
  {
    params::finalTime = 2.;
    rhoLeft       = 1.;         rhoRight      = 1.;
    pressureLeft  = 1.;         pressureRight = 1.;
    u1Left        = 0.;         u1Right       = 3.7;
    u2Left        = 0.;         u2Right       = 5.76;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 3.;         B1Right       = 3.;
    B2Left        = 3.;         B2Right       = -6.857;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "shock_tube_1")
  {
    params::finalTime = 1.;
    rhoLeft       = 1.;         rhoRight      = .1;
    pressureLeft  = 1000.;      pressureRight = 1.;
    u1Left        = 0.;         u1Right       = 0.;
    u2Left        = 0.;         u2Right       = 0.;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 1.;         B1Right       = 1.;
    B2Left        = 0.;         B2Right       = 0.;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "shock_tube_2")
  {
    params::finalTime = 1.;
    rhoLeft       = 1.;         rhoRight      = .1;
    pressureLeft  = 30.;        pressureRight = 1.;
    u1Left        = 0.;         u1Right       = 0.;
    u2Left        = 0.;         u2Right       = 0.;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 0.;         B1Right       = 0.;
    B2Left        = 20.;        B2Right       = 0.;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "collision")
  {
    params::finalTime      = 1.22;
    params::maxDtIncrement = 1.01;
    rhoLeft       = 1.;         rhoRight      = 1.;
    pressureLeft  = 1.;         pressureRight = 1.;
    u1Left        = 5.;         u1Right       = -5.;
    u2Left        = 0.;         u2Right       = 0.;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 10.;        B1Right       = 10.;
    B2Left        = 10.;        B2Right       = -10.;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (   params::shockTest == "stationary_shock" 
           || params::shockTest == "stationary_shock_BVP_input"
          )
  {
    rhoLeft       = 1.;         rhoRight      = 3.08312999;
    pressureLeft  = 1*(params::adiabaticIndex-1.); 
    pressureRight = 4.94577705*(params::adiabaticIndex-1.);
    u1Left        = 1.;         u1Right       = 0.32434571;
    u2Left        = 0.;         u2Right       = 0.;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 0.00001;    B1Right       = 0.00001;
    B2Left        = 0.;         B2Right       = 0.;
    B3Left        = 0.;         B3Right       = 0.;
  }

  af::seq domainLeftHalfX1(0, N1/2);
  af::seq domainRightHalfX1(N1/2, af::end);

  primOld->vars[vars::RHO](domainLeftHalfX1, span, span) = rhoLeft;
  primOld->vars[vars::U](domainLeftHalfX1, span, span)   = pressureLeft
                                                          /(params::adiabaticIndex-1.);
  primOld->vars[vars::U1](domainLeftHalfX1, span, span)  = u1Left;
  primOld->vars[vars::U2](domainLeftHalfX1, span, span)  = u2Left;
  primOld->vars[vars::U3](domainLeftHalfX1, span, span)  = u3Left;
  primOld->vars[vars::B1](domainLeftHalfX1, span, span)  = B1Left;
  primOld->vars[vars::B2](domainLeftHalfX1, span, span)  = B2Left;
  primOld->vars[vars::B3](domainLeftHalfX1, span, span)  = B3Left;

  primOld->vars[vars::RHO](domainRightHalfX1, span, span) = rhoRight;
  primOld->vars[vars::U](domainRightHalfX1, span, span)   = pressureRight
                                                          /(params::adiabaticIndex-1.);
  primOld->vars[vars::U1](domainRightHalfX1, span, span)  = u1Right;
  primOld->vars[vars::U2](domainRightHalfX1, span, span)  = u2Right;
  primOld->vars[vars::U3](domainRightHalfX1, span, span)  = u3Right;
  primOld->vars[vars::B1](domainRightHalfX1, span, span)  = B1Right;
  primOld->vars[vars::B2](domainRightHalfX1, span, span)  = B2Right;
  primOld->vars[vars::B3](domainRightHalfX1, span, span)  = B3Right;

  for (int var=0; var<vars::dof; var++)
  {
    primIC->vars[var] = primOld->vars[var];
  }

  if (params::shockTest == "stationary_shock_BVP_input")
  {
    PetscPrintf(PETSC_COMM_WORLD, "Setting up initial conditions...");

    int numGhost = params::numGhost;
    double xCoords1D[N1];
    double rho1D[N1];
    double u1D[N1];
    double u1_1D[N1];
    double q1D[N1];
    double dP1D[N1];

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank); /* get current proc id */

    if (rank==0)
    {
      FILE *rhoFile;
      FILE *uFile;
      FILE *u1File;
      FILE *qFile;
      FILE *dPFile;
      FILE *xCoordsFile;
  
      char *rhoLine     = NULL;
      char *uLine       = NULL;
      char *u1Line      = NULL; 
      char *qLine       = NULL; 
      char *dPLine      = NULL; 
      char *xCoordsLine = NULL;

      size_t rhoLen     = 0; ssize_t rhoRead;
      size_t uLen       = 0; ssize_t uRead;
      size_t u1Len      = 0; ssize_t u1Read;
      size_t qLen       = 0; ssize_t qRead;
      size_t dPLen      = 0; ssize_t dPRead;
      size_t xCoordsLen = 0; ssize_t xCoordsRead;

      rhoFile     = fopen(problemParams::rhoInputFile.c_str(), "r");
      uFile       = fopen(problemParams::uInputFile.c_str()  , "r");
      u1File      = fopen(problemParams::u1InputFile.c_str() , "r");
      qFile       = fopen(problemParams::qInputFile.c_str()  , "r");
      dPFile      = fopen(problemParams::dPInputFile.c_str() , "r");
      xCoordsFile = fopen(problemParams::xCoordsInputFile.c_str(), "r");

      if (   rhoFile      == NULL 
          || uFile        == NULL
          || u1File       == NULL
          || qFile        == NULL
          || dPFile       == NULL
          || xCoordsFile  == NULL
        )
      {
        PetscPrintf(PETSC_COMM_WORLD, "Input data files not found!\n");
        exit(1);
      }

      for (int i=0; i<N1; i++)
      {
        rhoRead     = getline(&rhoLine    , &rhoLen     , rhoFile);
        uRead       = getline(&uLine      , &uLen       , uFile);
        u1Read      = getline(&u1Line     , &u1Len      , u1File);
        qRead       = getline(&qLine      , &qLen       , qFile);
        dPRead      = getline(&dPLine     , &dPLen      , dPFile);
        xCoordsRead = getline(&xCoordsLine, &xCoordsLen , xCoordsFile);

        if (   rhoRead      == -1
            || uRead        == -1
            || u1Read       == -1
            || qRead        == -1
            || dPRead       == -1
            || xCoordsRead  == -1
          )
        {
          PetscPrintf(PETSC_COMM_WORLD,
          "Found the input data files but error in reading them! Check number of grid zones in python script\n");
          exit(1);
        }

        rho1D[i]     = atof(rhoLine);
        u1D[i]       = atof(uLine);
        u1_1D[i]     = atof(u1Line);
        q1D[i]       = atof(qLine);
        dP1D[i]      = atof(dPLine);
        xCoords1D[i] = atof(xCoordsLine);
      }

      free(rhoLine);
      free(uLine);
      free(u1Line);
      free(qLine);
      free(dPLine);
      free(xCoordsLine);

      fclose(rhoFile);
      fclose(uFile);
      fclose(u1File);
      fclose(qFile);
      fclose(dPFile);
      fclose(xCoordsFile);
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "file read complete\n");

    /* Broadcast the data from rank 0 proc to all other procs */
    MPI_Bcast(&rho1D[0]     , N1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&u1D[0]       , N1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&u1_1D[0]     , N1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&q1D[0]       , N1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dP1D[0]      , N1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&xCoords1D[0] , N1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
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

    array &rhoAnalytic = primIC->vars[vars::RHO];
    array &uAnalytic   = primIC->vars[vars::U];
    array &u1Analytic  = primIC->vars[vars::U1];

    array xCoords[3], XCoords[3];
    geomCenter->getXCoords(XCoords);
    geomCenter->XCoordsToxCoords(XCoords,xCoords);
    array x1 = xCoords[directions::X1];

    for (int i=0; i<primOld->N1Total - 2*numGhost; i++)
    {
      int iGlobal = i + primOld->iLocalStart;

      if (fabs(x1(i+numGhost, 0, 0).scalar<double>() - xCoords1D[iGlobal]) > 1e-10)
      {
        PetscPrintf(PETSC_COMM_WORLD,
                   "Mismatch in xCoords! Check x1 coords in python script. x1 grim = %.18f, x1 pBVP python solver = %.18f, i = %d\n",
                 x1(i+numGhost, 0, 0).scalar<double>(), xCoords1D[iGlobal], i
                 );
        exit(1);
      }
      rho(i+numGhost, span, span) = rho1D[iGlobal];
      u(i+numGhost, span, span)   = u1D[iGlobal];
      u1(i+numGhost, span, span)  = u1_1D[iGlobal];

      rhoAnalytic(i+numGhost, span, span) = rho1D[iGlobal];
      uAnalytic(i+numGhost, span, span)   = u1D[iGlobal];
      u1Analytic(i+numGhost, span, span)  = u1_1D[iGlobal];

      double pressure    = (params::adiabaticIndex - 1.)*u1D[iGlobal];
      double temperature = pressure/rho1D[iGlobal];
      double soundSpeed  = sqrt(params::adiabaticIndex*pressure
                                /(rho1D[iGlobal] 
                                  + params::adiabaticIndex*u1D[iGlobal]
                                 )
                               );
      if (params::conduction)
      {
        array &qAnalytic   = primIC->vars[vars::Q];
        qAnalytic(i+numGhost, span, span) = q1D[iGlobal];

        if (params::highOrderTermsConduction)
        {
          /* The evolved quantity is qTilde, but the initial data from the BVP
           * solver is given for q.  Need to rescale. Taken from physics.cpp */
          array &qTilde = primOld->vars[vars::Q];
        
          qTilde(i+numGhost, span, span) = 
            q1D[iGlobal] / temperature 
          / sqrt(rho1D[iGlobal] * params::ConductionAlpha * soundSpeed * soundSpeed);
        }
        else
        {
          array &q = primOld->vars[vars::Q];
          q(i+numGhost, span, span) = q1D[iGlobal];
        }
      }

      if (params::viscosity)
      {
        array &dPAnalytic  = primIC->vars[vars::DP];
        dPAnalytic(i+numGhost, span, span) = dP1D[iGlobal];

        if (params::highOrderTermsViscosity)
        {
          array &deltaPTilde = primOld->vars[vars::DP];

          deltaPTilde(i+numGhost, span, span) = 
            dP1D[iGlobal]
          / sqrt(  temperature *rho1D[iGlobal] * params::ViscosityAlpha 
                 * soundSpeed * soundSpeed
                );
        }
        else
        {
          array &dP = primOld->vars[vars::DP];
          dP(i+numGhost, span, span) = dP1D[iGlobal];
        }
      }
      primOld->vars[vars::B1] = 1e-5;
      
    }
  } /* End of "stationary_shock_BVP_input" */

  for (int var=0; var < numVars; var++)
  {
    primOld->vars[var].eval();
  }
  af::sync();

  fullStepDiagnostics(numReads, numWrites);
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{

}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
  bool WriteData  = 
    (floor(time/params::WriteDataEveryDt) != 
     floor((time-dt)/params::WriteDataEveryDt)
    );
  if(WriteData)
    {
      long long int WriteIdx = floor(time/params::WriteDataEveryDt);
      if(WriteIdx==0)
      {
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        PetscPrintf(PETSC_COMM_WORLD, "  Printing metric at zone CENTER...");
        geomCenter->gCovGrid->dump("gCov","gCov.h5");
        geomCenter->gConGrid->dump("gCon","gCon.h5");
        geomCenter->gGrid->dump("sqrtDetg","sqrtDetg.h5");
        geomCenter->xCoordsGrid->dump("xCoords","xCoords.h5");
        PetscPrintf(PETSC_COMM_WORLD, "done\n\n");
      }
      std::string filename = "primVarsT";
      std::string s_idx = std::to_string(WriteIdx);
      for(int i=0;i<6-s_idx.size();i++)
      {
        filename=filename+"0";
      }
      filename=filename+s_idx;
      filename=filename+".h5";

      primOld->dump("primitives",filename);
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
    }

  if (params::shockTest == "stationary_shock_BVP_input")
  {
    elemOld->set(*primOld, *geomCenter, numReads, numWrites);
    double errorRho = 
      af::norm(af::flat(elemOld->rho - primIC->vars[vars::RHO])(domainX1, span, span));

    double errorU   = 
      af::norm(af::flat(elemOld->u   - primIC->vars[vars::U])(domainX1, span, span));

    double errorU1  = 
      af::norm(af::flat(elemOld->u1  - primIC->vars[vars::U1])(domainX1, span, span));

    double errorQ   = 
      af::norm(af::flat(elemOld->q   - primIC->vars[vars::Q])(domainX1, span, span));

    double errordP  = 
      af::norm(af::flat(elemOld->deltaP  - primIC->vars[vars::DP])(domainX1, span, span));

    errorRho = errorRho/N1/N2/N3;
    errorU   = errorU/N1/N2/N3;
    errorU1  = errorU1/N1/N2/N3;
    errorQ   = errorQ/N1/N2/N3;
    errordP  = errordP/N1/N2/N3;

    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "    ---EMHD Shock Diagnostics---\n");
    PetscPrintf(PETSC_COMM_WORLD, "     Error in rho = %e\n", errorRho);
    PetscPrintf(PETSC_COMM_WORLD, "     Error in u   = %e\n", errorU  );
    PetscPrintf(PETSC_COMM_WORLD, "     Error in u1  = %e\n", errorU1 );
    PetscPrintf(PETSC_COMM_WORLD, "     Error in q   = %e\n", errorQ  );
    PetscPrintf(PETSC_COMM_WORLD, "     Error in dP  = %e\n", errordP );
  }
}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{
  if (params::shockTest == "stationary_shock_BVP_input")
  {
    int numGhost = params::numGhost;
    // 1) Choose which primitive variables are corrected.
    grid* primBC;
    if(currentStep == timeStepperSwitches::HALF_STEP)
    {
      primBC = primOld;
    }
    else
    {
      primBC = primHalfStep;
    }
    af::seq leftBoundary(0, numGhost-2);
    af::seq rightBoundary(primOld->N1Local + numGhost+2, 
                          primOld->N1Local + 2*numGhost-1
                         );

    if (primOld->iLocalStart == 0)
    {
      for (int var=0; var < vars::dof; var++)
      {
        primBC->vars[var](leftBoundary, span, span) = 
          primIC->vars[var](leftBoundary, span, span);
      }
    }

    if (primOld->iLocalStart == N1)
    {
      for (int var=0; var < vars::dof; var++)
      {
        primBC->vars[var](rightBoundary, span, span) = 
          primIC->vars[var](rightBoundary, span, span);
      }
    }
  }

}

void timeStepper::applyProblemSpecificFluxFilter(int &numReads,int &numWrites)
{

}

int timeStepper::CheckWallClockTermination()
{
  return 0;
}
