add_library(grid grid.cpp grid.hpp)

set_source_files_properties(gridPy.pyx PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module(gridPy gridPy.pyx grid.cpp ../params.cpp)
target_link_libraries(gridPy ${PETSC_LIBRARIES} ${ArrayFire_CPU_LIBRARIES})
