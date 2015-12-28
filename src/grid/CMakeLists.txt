add_library(grid grid.cpp grid.hpp)

#set_source_files_properties(gridpy.pyx PROPERTIES CYTHON_IS_CXX TRUE)

#cython_add_module(gridpy gridpy.pyx grid.cpp)
#target_link_libraries(gridpy ${PETSC_LIBRARIES} ${ArrayFire_CPU_LIBRARIES})
