#include "grim.hpp"
#include "params.cpp"

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
    timeStepper ts;

    /*af::dtype type;
    array arr = (ts.geom->xCoords[locations::CENTER][1]);
    type = arr.type();
    printf("xCoords has type %i\n",type);*/
    

    double Aw = 1.e-5;
    double k1 = 2.*M_PI;
    double k2 = 4.*M_PI;
    array cphi = af::cos(k1*ts.geom->xCoords[locations::CENTER][1]
			 +k2*ts.geom->xCoords[locations::CENTER][2]);
    array sphi = af::sin(k1*ts.geom->xCoords[locations::CENTER][1]
			 +k2*ts.geom->xCoords[locations::CENTER][2]);

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
    ts.primOld->vars[vars::RHO] = 1.;
    ts.primOld->vars[vars::U]   = 2.;
    ts.primOld->vars[vars::U1]  = 0.;
    ts.primOld->vars[vars::U2]  = 0.; 
    ts.primOld->vars[vars::U3]  = 0.; 
    ts.primOld->vars[vars::B1]  = 0.1; 
    ts.primOld->vars[vars::B2]  = 0.3;
    ts.primOld->vars[vars::B3]  = 0.;
    ts.primOld->vars[vars::Q]  = 0.;
    ts.primOld->vars[vars::DP]  = 0.;
    ts.primOld->vars[vars::RHO] += Aw*cphi*(-0.518522524082246)
      +Aw*sphi*0.1792647678001878;
    ts.primOld->vars[vars::U] += Aw*cphi*0.5516170736393813;
    ts.primOld->vars[vars::U1]+= Aw*cphi*0.008463122479547856
      +Aw*sphi*(-0.011862022608466367);
    ts.primOld->vars[vars::U2]+= Aw*cphi*(-0.16175466371870734)
      +Aw*sphi*(0.034828080823603294);
    ts.primOld->vars[vars::B1]+= Aw*cphi*(-0.05973794979640743)
      +Aw*sphi*0.03351707506150924;
    ts.primOld->vars[vars::B2]+=Aw*cphi*0.02986897489820372
      -Aw*sphi*0.016758537530754618;
    ts.primOld->vars[vars::Q]+=Aw*cphi*0.5233486841539436
      -Aw*sphi*0.04767672501939603;
    ts.primOld->vars[vars::DP]+=Aw*cphi*0.2909106062057657
      -Aw*sphi*0.02159452055336572;
    double Gamma = - 0.5533585207638141;
    double Omega = - 3.6262571286888425;

    params::Time = 0.;
    while(params::Time<0.5)
      {
	ts.timeStep();
	params::Time+=params::dt;
	
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
	cphi = af::cos(k1*ts.geom->xCoords[locations::CENTER][1]
		       +k2*ts.geom->xCoords[locations::CENTER][2]
		       +Omega*params::Time);
	sphi = af::sin(k1*ts.geom->xCoords[locations::CENTER][1]
		       +k2*ts.geom->xCoords[locations::CENTER][2]
		       +Omega*params::Time);
	array rhoan = 1.+(Aw*cphi*(-0.518522524082246)
			  +Aw*sphi*0.1792647678001878)*exp(Gamma*params::Time);
	array psian = 0. + (Aw*cphi*0.2909106062057657
			    -Aw*sphi*0.02159452055336572)*exp(Gamma*params::Time);
	double error = af::norm(af::flat((ts.primOld->vars[vars::RHO]
	- rhoan)));
	double error2 = af::norm(af::flat((ts.primOld->vars[vars::DP]
					   -psian)));

	error = error/params::N1/params::N2/params::N3;
	error2 = error2/params::N1/params::N2/params::N3;
	printf("Time = %e; dt = %e; Error = %e; Error2 = %e\n",params::Time,params::dt,error,error2);
	//af_print(ts.primOld->vars[vars::RHO](span,0,0)-rhoan(span,0,0),12);
      }

  }

  PetscFinalize();  
  return(0);
}
