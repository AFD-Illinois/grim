#ifndef _GRIM_H
#define _GRIM_H

//#define VIENNACL_BUILD_INFO
//#define VIENNACL_WITH_OPENCL

#include <petsc.h>
#include <petscdmda.h>
//#include <petscviennacl.h>
#include <petscviewerhdf5.h>
//#include <CL/cl.hpp>
#include <fstream>
#include <ctime>
#include "constants.h"

static const char help[] = 
    "GRIM -- General Relativistic Implicit Magnetohydrodynamics";

//cl_int clErr;
//std::vector<cl::Platform> platforms;
//std::vector<cl::Device> devices;
//cl::Context context;
//cl::CommandQueue queue;
//cl::Program program;
//cl::Kernel kernel;

//viennacl::ocl::program program;
//
//void CheckCLErrors(cl_int clErr, std::string errMsg)
//{
//    if (clErr!=CL_SUCCESS) {
//        SETERRQ(PETSC_COMM_SELF, 1, errMsg.c_str());
//    }
//}
extern void InitialCondition(TS ts, Vec prim);
extern void Benchmark(TS ts, Vec prim);
extern PetscErrorCode ComputeResidual(TS ts,
                               PetscScalar t,
                               Vec Prim, Vec dPrim_dt,
                               Vec F, void *ptr);

#endif
