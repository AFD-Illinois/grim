add_library(boundary boundary.cpp boundary.hpp)
target_link_libraries(boundary grid)

#set_source_files_properties(boundaryPy.pyx
#                            PROPERTIES CYTHON_IS_CXX TRUE)

#cython_add_module(boundaryPy boundaryPy.pyx)

#target_link_libraries(boundaryPy boundary
#                      ${ArrayFire_LIBRARIES}
#                     )
