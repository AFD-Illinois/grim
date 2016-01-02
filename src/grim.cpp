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
    

    double Aw = 1.e-4;
    array cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]);
    array sphi = af::sin(2*M_PI*ts.geom->xCoords[locations::CENTER][1]);

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
    ts.primOld->vars[vars::RHO] = 1.+Aw*0.408365507885*cphi;
    ts.primOld->vars[vars::U]   = 2.+Aw*0.816299597519*cphi;
    ts.primOld->vars[vars::U1]  = 0.+Aw*0.0054163532418*sphi; 
    ts.primOld->vars[vars::U2]  = 0.; 
    ts.primOld->vars[vars::U3]  = 0.; 
    ts.primOld->vars[vars::B1]  = 0.01; 
    ts.primOld->vars[vars::B2]  = 0.;
    ts.primOld->vars[vars::B3]  = 0.;
    ts.primOld->vars[vars::Q]  = 0.-Aw*0.00361662427435*sphi;
    ts.primOld->vars[vars::DP]  = 0.+Aw*0.408472963863*cphi;

    params::Time = 0.;
    while(params::Time<0.05)
      {
	ts.timeStep(params::dt);
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
	cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
		       -0.0833369872094*params::Time);
	sphi = af::sin(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
		       -0.0833369872094*params::Time);
	array rhoan = 1.+Aw*0.408365507885*cphi;
	double error = af::norm(af::flat((ts.primOld->vars[vars::RHO]
	- rhoan)));

	error = error/params::N1/params::N2/params::N3;
	printf("Time = %e; dt = %e; Error = %e\n",params::Time,params::dt,error);
	//af_print(ts.primOld->vars[vars::RHO](span,0,0)-rhoan(span,0,0),12);
      }

  }

  PetscFinalize();  
  return(0);
}
