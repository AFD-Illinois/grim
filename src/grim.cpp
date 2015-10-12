#include "grim.h"
#include <petscviennacl.h>
#include <boost/preprocessor/stringize.hpp>
//#include "inputs.h"
#define VEX_STRINGIZE_SOURCE(...) #__VA_ARGS__

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_size(PETSC_COMM_WORLD, &rank);

  typedef std::vector< viennacl::ocl::platform > platforms_type;
  platforms_type platforms = viennacl::ocl::get_platforms();
  bool is_first_platform = true;
  for (platforms_type::iterator platform_iter  = platforms.begin();
                                platform_iter != platforms.end();
                              ++platform_iter
      )
  {
    typedef std::vector<viennacl::ocl::device> devices_type;
    devices_type devices = platform_iter->devices(CL_DEVICE_TYPE_ALL);
    FILE *opencl_info; 
    opencl_info = fopen("opencl_info.txt", "w");
    PetscSynchronizedFPrintf
      (PETSC_COMM_WORLD, opencl_info,
       "=============================================\n"
       "  Platform Information for MPI proc %d of %d \n"
       "=============================================\n"
       "                                             \n"
       "Vendor and version: %s                       \n"
       ,rank, size, platform_iter->info().c_str()
      );

    if (is_first_platform)
    {
      PetscSynchronizedFPrintf
        (PETSC_COMM_WORLD, opencl_info,
         "Using the above OpenCL platform by default.\n\n"
        );
      is_first_platform = false;
    }

    PetscSynchronizedFPrintf
      (PETSC_COMM_WORLD, opencl_info,
       "Available devices:\n\n"
      );

    int device_num=1;
    for (devices_type::iterator iter = devices.begin();
                                iter != devices.end(); 
                                iter++
        )
    {
      PetscSynchronizedFPrintf
        (PETSC_COMM_WORLD, opencl_info,
         "=============================================\n"
         "Device %d\n\n"
         "%s         \n"
         "=============================================\n",
         device_num, iter->full_info().c_str()
        );
      device_num++;
    }
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD, opencl_info);
  }
  
  //struct timeStepper ts;
  //timeStepperInit(&ts);

//  while (ts.t + ts.dt < FINAL_TIME)
//  {
//    timeStep(&ts);
//  }
//  
//  /* One final step */
//  if (ts.t < FINAL_TIME)
//  {
//    ts.dt = FINAL_TIME - ts.t;
//    timeStep(&ts);
//  }
//
//  timeStepperDestroy(&ts);


  DM dm;

  DMBoundaryType boundaryX1 = DM_BOUNDARY_GHOSTED;
  DMBoundaryType boundaryX2 = DM_BOUNDARY_GHOSTED;

  int numX1 = 128;
  int numX2 = 128;
  int numVar = 10;

  DMDACreate2d(PETSC_COMM_WORLD, 
               boundaryX1, boundaryX2,
               DMDA_STENCIL_BOX,
               numX1, numX2,
               PETSC_DECIDE, PETSC_DECIDE,
               numVar, 0,
               PETSC_NULL, PETSC_NULL,
               &dm
              );


  PetscFinalize();  
  return(0);
}
