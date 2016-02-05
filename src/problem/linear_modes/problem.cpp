#include "../problem.hpp"

double Aw = 1.e-5;
double k1 = 2.*M_PI;
double k2 = 4.*M_PI;
double Gamma = - 0.5533585207638141;
double Omega = - 3.6262571286888425;

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
}

void timeStepper::initialConditions(const array xCoords[3],
                                    array prim[vars::dof]
                                   )
{

  array cphi = af::cos(  k1*xCoords[directions::X1]
                  		 + k2*xCoords[directions::X2]
                      );

  array sphi = af::sin(  k1*xCoords[directions::X1]
                  		 + k2*xCoords[directions::X2]
                      );

  /* Initial conditions */

  //Alfven Wave
  /*prim[vars::RHO] = 1.;
  prim[vars::U]   = 2.;
  prim[vars::U1]  = 0.;
  prim[vars::U2]  = Aw*0.462905090215*cphi;
  prim[vars::U3]  = 0.;
  prim[vars::B1]  = 0.01;
  prim[vars::B2]  = Aw*0.886407850514*cphi;
  prim[vars::B3]  = 0.;*/

  //Sound wave
  /*prim[vars::RHO] = 1.+Aw*0.345991032308*cphi;
  prim[vars::U]   = 2.+Aw*0.922642752822*cphi;
  prim[vars::U1]  = 0.-Aw*0.170354208129*cphi; 
  prim[vars::U2]  = 0.; 
  prim[vars::U3]  = 0.; 
  prim[vars::B1]  = 0.01; 
  prim[vars::B2]  = 0.;
  prim[vars::B3]  = 0.;*/

  //Full EMHD mode (from grim2D)
  prim[vars::RHO] = 1.;
  prim[vars::U]   = 2.;
  prim[vars::U1]  = 0.;
  prim[vars::U2]  = 0.; 
  prim[vars::U3]  = 0.; 
  prim[vars::B1]  = 0.1; 
  prim[vars::B2]  = 0.3;
  prim[vars::B3]  = 0.;
  prim[vars::Q]   = 0.;
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

  for (int var=0; var < numVars; var++)
  {
    prim[var].eval();
  }
  af::sync();
}

void timeStepper::halfStepDiagnostics(const array xCoords[3],
                                      array prim[vars::dof]
                                     )
{

}

void timeStepper::fullStepDiagnostics(const array xCoords[3],
                                      array prim[vars::dof]
                                     )
{
  /* Compute the errors for the different modes */

  //Alfven wave
	/*cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
		       +0.0328124176673*params::Time);
	array u2an = Aw*0.462905090215*cphi;
	double error = af::norm(af::flat((ts.primOld->vars[vars::U2]
	- u2an)));*/

	//MHD Sound wave
	/*cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
	  +3.09362659024*params::Time);
	array rhoan = 1.+Aw*0.345991032308*cphi;
	double error = af::norm(af::flat((ts.primOld->vars[vars::RHO]
	- rhoan)));*/

	//EMHD Sound wave
	array cphi = af::cos(  k1*xCoords[directions::X1]
      	      	       + k2*xCoords[directions::X2]
		                   + Omega*time
                      );

	array sphi = af::sin(  k1*xCoords[directions::X1]
		                   + k2*xCoords[directions::X2]
		                   + Omega*time
                      );

	array rhoAnalytic = 1. + (  Aw*cphi*(-0.518522524082246)
			                      + Aw*sphi*0.1792647678001878
                           )*exp(Gamma*time);

	array dPAnalytic = 0. + ( Aw*cphi*0.2909106062057657
			                     -Aw*sphi*0.02159452055336572
                          )*exp(Gamma*time);

	double errorRho = af::norm(af::flat(prim[vars::RHO] - rhoAnalytic));

	double errordP  = af::norm(af::flat(prim[vars::DP]  - dPAnalytic));

	//af_print(prim[vars::RHO] - rhoAnalytic,12);

	errorRho = errorRho/N1/N2/N3;
	errordP  = errordP/N1/N2/N3;
	PetscPrintf(PETSC_COMM_WORLD, 
              "Time = %e; dt = %e; Error in rho = %e; Error in dP = %e\n",
              time, dt, errorRho, errordP
              );
}

void timeStepper::setProblemSpecificBCs(const array xCoords[3],
                                        array prim[vars::dof]
                                       )
{

};
