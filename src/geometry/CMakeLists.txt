add_library(geometry geometry.cpp geometry.hpp)
target_link_libraries(geometry grid)

set_source_files_properties(geometryPy.pyx PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module(geometryPy geometryPy.pyx)
target_link_libraries(geometryPy geometry ${ArrayFire_LIBRARIES})
