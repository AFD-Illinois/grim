#include "grim.h"
#include <petscviennacl.h>
#include <CL/cl.hpp>
#include <yaml-cpp/yaml.h>

#include <boost/preprocessor/stringize.hpp>

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

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


  if (rank==0)
  {
    YAML::Node input = YAML::LoadFile("input.yaml");

    printf("N1 = %d, N2 = %d\n", 
           input["domain"]["N1"].as<int>(),
           input["domain"]["N2"].as<int>()
          );

    printf("left = %s, right = %s\n",
           input["domain"]["boundaries"]["left"].as<std::string>().c_str(),
           input["domain"]["boundaries"]["right"].as<std::string>().c_str()
          );

  }

static const char *my_compute_program =
BOOST_PP_STRINGIZE(
  __kernel void elementwise_prod(__global const float * vec1,
                                 __global const float * vec2,
                                 __global float * result,
                                 unsigned int size
                                )
  {
    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))
      result[i] = vec1[i] * vec2[i];
  };
);

   viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(my_compute_program, "my_compute_program");

  viennacl::ocl::kernel & my_kernel_mul = my_prog.get_kernel("elementwise_prod");

  viennacl::vector<double> vec1(numX1);
  viennacl::vector<double> vec2(numX1);
  viennacl::vector<double> result_mul(numX1);

  cl::Context context;
   viennacl::ocl::enqueue(my_kernel_mul(vec1, vec2, result_mul, static_cast<cl_uint>(vec1.size())));

  PetscFinalize();  
  return(0);
}
