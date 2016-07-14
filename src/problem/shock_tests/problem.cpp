#include "../problem.hpp"

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
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
    params::finalTime = 1.22;
    rhoLeft       = 1.;         rhoRight      = 1.;
    pressureLeft  = 1.;         pressureRight = 1.;
    u1Left        = 5.;         u1Right       = -5.;
    u2Left        = 0.;         u2Right       = 0.;
    u3Left        = 0.;         u3Right       = 0.;
    B1Left        = 10.;        B1Right       = 10.;
    B2Left        = 10.;        B2Right       = -10.;
    B3Left        = 0.;         B3Right       = 0.;
  }
  else if (params::shockTest == "stationary_shock")
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

}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{

}

void timeStepper::applyProblemSpecificFluxFilter(int &numReads,int &numWrites)
{

}
