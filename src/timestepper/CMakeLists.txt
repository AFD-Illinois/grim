add_library(timestepper timestepper.cpp timestepper.hpp timestep.cpp 
            fvmfluxes.cpp residual.cpp solve.cpp constrainedtransport.cpp)
target_link_libraries(timestepper geometry grid physics)

#set_source_files_properties(timeStepperPy.pyx PROPERTIES CYTHON_IS_CXX TRUE)
#cython_add_module(timeStepperPy timeStepperPy.pyx)
#target_link_libraries(timeStepperPy timestepper problem
#                      reconstruction boundary ${PETSC_LIBRARIES}
#                      )
