add_library(reconstruction reconstruction.hpp minmod.cpp weno5.cpp
            reconstruction.cpp)

set_source_files_properties(reconstructionPy.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module(reconstructionPy reconstructionPy.pyx reconstruction.cpp
                  weno5.cpp minmod.cpp ../params.cpp
                 )
target_link_libraries(reconstructionPy
                      ${ArrayFire_CPU_LIBRARIES}
                     )
