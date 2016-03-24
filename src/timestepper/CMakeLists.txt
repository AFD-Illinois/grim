add_library(timestepper timestepper.cpp timestepper.hpp timestep.cpp 
            fvmfluxes.cpp residual.cpp solve.cpp constrainedtransport.cpp)
target_link_libraries(timestepper geometry grid physics)

set_source_files_properties(timestepperPy.pyx PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(timeStepperPy timestepperPy.pyx
                  ${CMAKE_SOURCE_DIR}/problem/${PROBLEM}/params.cpp
                 )
target_link_libraries(timeStepperPy timestepper reconstruction ${ArrayFire_LIBRARIES})
