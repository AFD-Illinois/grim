add_library(reconstruction reconstruction.hpp minmod.cpp weno5.cpp ppm.cpp
            reconstruction.cpp)
target_link_libraries(reconstruction ${ArrayFire_LIBRARIES})

set_source_files_properties(reconstructionPy.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module(reconstructionPy
                  reconstructionPy.pyx
                  )

target_link_libraries(reconstructionPy reconstruction params)
