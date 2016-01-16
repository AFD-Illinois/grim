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

void timeStepper::initialConditions()
{
  array cphi = af::cos(  k1*geom->xCoords[locations::CENTER][1]
                  		 + k2*geom->xCoords[locations::CENTER][2]
                      );

  array sphi = af::sin(  k1*geom->xCoords[locations::CENTER][1]
	                  	 + k2*geom->xCoords[locations::CENTER][2]
                      );

  /* Initial conditions */

  //Alfven Wave
  /*ts.primOld->vars[vars::RHO] = 1.;
  ts.primOld->vars[vars::U]   = 2.;
  ts.primOld->vars[vars::U1]  = 0.;
  ts.primOld->vars[vars::U2]  = Aw*0.462905090215*cphi;
  ts.primOld->vars[vars::U3]  = 0.;
  ts.primOld->vars[vars::B1]  = 0.01;
  ts.primOld->vars[vars::B2]  = Aw*0.886407850514*cphi;
  ts.primOld->vars[vars::B3]  = 0.;*/

  //Sound wave
  /*ts.primOld->vars[vars::RHO] = 1.+Aw*0.345991032308*cphi;
  ts.primOld->vars[vars::U]   = 2.+Aw*0.922642752822*cphi;
  ts.primOld->vars[vars::U1]  = 0.-Aw*0.170354208129*cphi; 
  ts.primOld->vars[vars::U2]  = 0.; 
  ts.primOld->vars[vars::U3]  = 0.; 
  ts.primOld->vars[vars::B1]  = 0.01; 
  ts.primOld->vars[vars::B2]  = 0.;
  ts.primOld->vars[vars::B3]  = 0.;*/

  //Full EMHD mode (from grim2D)
  primOld->vars[vars::RHO] = 1.;
  primOld->vars[vars::U]   = 2.;
  primOld->vars[vars::U1]  = 0.;
  primOld->vars[vars::U2]  = 0.; 
  primOld->vars[vars::U3]  = 0.; 
  primOld->vars[vars::B1]  = 0.1; 
  primOld->vars[vars::B2]  = 0.3;
  primOld->vars[vars::B3]  = 0.;
  primOld->vars[vars::Q]  = 0.;
  primOld->vars[vars::DP]  = 0.;
  primOld->vars[vars::RHO] += Aw*cphi*(-0.518522524082246)
    +Aw*sphi*0.1792647678001878;
  primOld->vars[vars::U] += Aw*cphi*0.5516170736393813;
  primOld->vars[vars::U1]+= Aw*cphi*0.008463122479547856
    +Aw*sphi*(-0.011862022608466367);
  primOld->vars[vars::U2]+= Aw*cphi*(-0.16175466371870734)
    +Aw*sphi*(0.034828080823603294);
  primOld->vars[vars::B1]+= Aw*cphi*(-0.05973794979640743)
    +Aw*sphi*0.03351707506150924;
  primOld->vars[vars::B2]+=Aw*cphi*0.02986897489820372
    -Aw*sphi*0.016758537530754618;
  primOld->vars[vars::Q]+=Aw*cphi*0.5233486841539436
    -Aw*sphi*0.04767672501939603;
  primOld->vars[vars::DP]+=Aw*cphi*0.2909106062057657
    -Aw*sphi*0.02159452055336572;

  for (int var=0; var<vars::dof; var++)
  {
    primOld->vars[var].eval();
  }
  af::sync();
}

void timeStepper::halfStepDiagnostics()
{

}

void timeStepper::fullStepDiagnostics()
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
	array cphi = af::cos(  k1*geom->xCoords[locations::CENTER][1]
      	      	       + k2*geom->xCoords[locations::CENTER][2]
		                   + Omega*params::Time
                      );

	array sphi = af::sin(  k1*geom->xCoords[locations::CENTER][1]
		                   + k2*geom->xCoords[locations::CENTER][2]
		                   + Omega*params::Time
                      );

	array rhoan = 1. + (  Aw*cphi*(-0.518522524082246)
			                + Aw*sphi*0.1792647678001878
                     )*exp(Gamma*params::Time);

	array psian = 0. + ( Aw*cphi*0.2909106062057657
			                -Aw*sphi*0.02159452055336572
                     )*exp(Gamma*params::Time);

	double error = af::norm( af::flat((  primOld->vars[vars::RHO]
	                                   - rhoan
                                    )
                                   )
                         );

	double error2 = af::norm(af::flat((  primOld->vars[vars::DP]
		                        			   - psian
                                    )
                                   )
                          );

	error = error/params::N1/params::N2/params::N3;
	error2 = error2/params::N1/params::N2/params::N3;
	printf("Time = %e; dt = %e; Error = %e; Error2 = %e\n",
         params::Time,params::dt,error,error2
        );
	//af_print(ts.primOld->vars[vars::RHO](span,0,0)-rhoan(span,0,0),12);
}

void timeStepper::setProblemSpecificBCs()
{

};
