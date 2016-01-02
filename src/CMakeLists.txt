# To compile the code, make a seperate build folder and then do the following:
#
# cd /location/to/build
# cmake /location/to/grim/src
# make
#
# If a different compiler than the one cmake detects is needed, then the
# recommended way is to do the following:
#
# cmake -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpic++ /location/to/grim/src 

cmake_minimum_required(VERSION 2.8)

project(grim)

# ------------------------------USER OPTIONS----------------------------------#
# Options for icc:
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3 -vec-threshold0 -qopt-report=5 -qopt-report-phase=vec -xhost")
# Options for gcc:
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x")

# Set custom install folders here
set(PETSC_DIR "/global/home/users/foucart/petsc_installation")
#set(YAML_ROOT_DIR "/global/home/users/manic/yaml-cpp/build")
set(ArrayFire_DIR "/global/home/users/foucart/arrayfire-no-gl-3.2.1/share/ArrayFire/cmake")

# Double or float precision?
set(REAL "double")

# Choose the problem to be solved:
#
# 1) linear_modes
# 2) shock_tests
# 3) advection_test
# 4) magnetized_field_loop_advection
# 5) orzag_tang
# 6) magnetized_explosion
# 7) atmosphere
# 8) bondi_inflow
# 10) anisotropic_conduction
# 11) torus
# 12) buoyancy_instabilities
# 13) bondi_viscous

set(PROBLEM "performance_test")
# Now choose problem specific options in the respective problem folder in
# src/problem/
# ------------------------------END OF USER OPTIONS---------------------------#

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Need to find all packages before starting compilation, so that all headers
# needed by the code are included.

# Should the external libraries needed by the internal libraries go into their
# own CMakeLists in their own folders?

# TODO: Only find packages depending on what options are set.

# System libraries
find_library(MATH_LIBRARIES m REQUIRED)

# External packages
find_package(PETSc REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)
include(UseCython)
#find_package(YAML REQUIRED)
find_package(ArrayFire REQUIRED)
if (ArrayFire_CPU_FOUND)
  message("ArrayFire_CPU_LIBRARIES: " ${ArrayFire_CPU_LIBRARIES})
endif()
if (ArrayFire_CUDA_FOUND)
  message("ArrayFire_CUDA_LIBRARIES: " ${ArrayFire_CUDA_LIBRARIES})
endif()
if (ArrayFire_OpenCL_FOUND)
  message("ArrayFire_OpenCL_LIBRARIES: " ${ArrayFire_OpenCL_LIBRARIES})
endif()

include_directories(${PETSC_INCLUDES} 
#                    ${YAML_INCLUDES}
                    ${ArrayFire_INCLUDE_DIRS}
                    ${PYTHON_NUMPY_INCLUDE_DIR}
                    grid
                   )

add_subdirectory(grid)
add_subdirectory(geometry)
add_subdirectory(physics)
add_subdirectory(reconstruction)
#add_subdirectory(boundary)
add_subdirectory(timestepper)

set(PROBLEM_DIR ${CMAKE_SOURCE_DIR}/problem/${PROBLEM})
include(${PROBLEM_DIR}/CMakeLists.txt)
#add_subdirectory(problem)

add_executable(grim grim.cpp grim.hpp params.hpp)

target_link_libraries(grim grid geometry reconstruction physics timestepper
                      ${MATH_LIBRARIES} 
                      ${PETSC_LIBRARIES}
#                      ${YAML_LIBRARIES}
                      ${ArrayFire_CPU_LIBRARIES}
                     )


message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("C Compiler       : " ${CMAKE_C_COMPILER})
message("C_FLAGS          : " ${CMAKE_C_FLAGS})
message("PETSc dir        : " ${PETSC_DIR})
message("")
