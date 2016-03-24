add_library(timestepper timestepper.cpp timestepper.hpp timestep.cpp 
            fvmfluxes.cpp residual.cpp solve.cpp constrainedtransport.cpp)
target_link_libraries(timestepper geometry grid physics)
