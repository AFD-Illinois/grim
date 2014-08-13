# To compiler the code, make a seperate build folder and then do the following:
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
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")

# Double or float precision?
set(REAL "double")

# Use as ARRAY_ARGS ARRAY_SIZE while declaring arrays in functions arguments
# Must use -std=c99 if using the "restrict" keyword while compiling 
set(ARRAY_ARGS "const restrict static")

# Set custom install folders here
set(PETSC_DIR "/home/mc/Downloads/petsc_optimized/")

# Dump folder
set(OUTPUT_DIR "/home/mc/PhD/grim_petsc/build/")

# Choose the problem to be solved here:
#
# 1) LINEAR_MODES_1D
# 2) KOMISSAROV_FAST_SHOCK_1D
# 3) KOMISSAROV_SLOW_SHOCK_1D
# 4) KOMISSAROV_SWITCH_OFF_1D
# 5) KOMISSAROV_SWITCH_ON_1D
# 6) KOMISSAROV_SHOCK_TUBE_1_1D
# 7) KOMISSAROV_SHOCK_TUBE_2_1D
# 8) KOMISSAROV_COLLISION_1D
# 9) MAGNETIZED_FIELD_LOOP_ADVECTION_2D
# 10) ORZAG_TANG_2D
# 11) KOMISSAROV_MAGNETIZED_CYLINDRICAL_EXPLOSION_2D
# 12) EQUILIBRIUM_TORUS_2D
# 13) TORUS_2D
# 14) MTI_TEST_2D
set(PROBLEM "MAGNETIZED_FIELD_LOOP_ADVECTION_2D")

# Time stepping options: EXPLICIT, IMEX or IMPLICIT
set(TIME_STEPPING "EXPLICIT")

# Adiabatic index for the equation of state
set(ADIABATIC_INDEX "4./3")

# Conduction: ON or OFF?
set(CONDUCTION "OFF")

# Domain size. If the problem is 1D then N2 is ignored.
set(N1 "64")
set(N2 "64")

# The entire global domain is divided into tiles that are small enough to fit
# into the cache of the compute node or an accelerator. This technique
# optimizes cache usage and makes explicit use of the deep memory hierarchies
# which are prevalent in all kinds of processors and accelerators. 
# 
# Caution: 1) The tile sizes need to divide the global domain size in each
#             direction exactly!
#          2) We want a tile size that is as large as the cache. Too large a
#             size and the code can fail silently because local memory in the
#             OpenCL specification is a different amount for different 
#             machines.
#             NEED TO PUT A CHECK FOR THIS.
set(TILE_SIZE_X1 "8")
set(TILE_SIZE_X2 "8")

# ------------------------------END OF USER OPTIONS---------------------------#

configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/inputs.h)

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
find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${PETSC_INCLUDES} ${GSL_INCLUDE_DIRS})

add_subdirectory(geometry)
add_subdirectory(physics)
add_subdirectory(reconstruct)
add_subdirectory(riemannsolver)
add_subdirectory(boundary)
add_subdirectory(gridzone)
add_subdirectory(timestepper)

add_executable(grim grim.c grim.h)

target_link_libraries(grim geometry physics reconstruct boundary gridzone
                      timestepper
                      ${MATH_LIBRARIES}
                      ${PETSC_LIBRARIES})

message("")
message("#########################")
message("# Configuration options #")
message("#########################")
message("")
message("C Compiler       : " ${CMAKE_C_COMPILER})
message("C_FLAGS          : " ${CMAKE_C_FLAGS})
message("PETSc dir        : " ${PETSC_DIR})
message("Output dir       : " ${OUTPUT_DIR})
message("Problem          : " ${PROBLEM})
message("Time stepping    : " ${TIME_STEPPING})
message("Adiabatic index  : " ${ADIABATIC_INDEX})
message("Heat conduction  : " ${CONDUCTION})
message("Resolution       : " ${N1} "x" ${N2})
message("Tile size        : " ${TILE_SIZE_X1} "x" ${TILE_SIZE_X2})
message("")
