# To compiler the code, make a seperate build folder and then do the following:
#
# cd /location/to/build
# cmake /location/to/reaper/src
# make
#
# If a different compiler than the one cmake detects is needed, then the
# recommended way is to do the following:
#
# cmake -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpic++ /location/to/reaper/src 

cmake_minimum_required(VERSION 2.8)

project(reaper)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")

# Set custom install folders here
set(PETSC_DIR "/home/mc/Downloads/petsc_optimized/")
set(CUBATURE_DIR "/home/mc/Downloads/cubature-1.0/")

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
find_package(GSL REQUIRED)
find_package(Cubature REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${PETSC_INCLUDES} ${GSL_INCLUDE_DIRS}
                    ${CUBATURE_INCLUDE_DIR})

add_subdirectory(geometry)
add_subdirectory(physics)
add_subdirectory(reconstruct)
add_subdirectory(boundary)
add_subdirectory(gridzone)
add_subdirectory(timestepper)

add_executable(reaper reaper.c reaper.h
               ${CUBATURE_INCLUDE_DIR}/hcubature.c)

target_link_libraries(reaper geometry physics reconstruct boundary gridzone
                      timestepper
                      ${MATH_LIBRARIES}
                      ${PETSC_LIBRARIES} ${GSL_LIBRARIES})
