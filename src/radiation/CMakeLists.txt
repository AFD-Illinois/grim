add_library(radiation radiation.hpp
            ${CMAKE_SOURCE_DIR}/problem/${PROBLEM}/params.cpp
            ${CMAKE_SOURCE_DIR}/problem/${PROBLEM}/problem.cpp
           )
target_link_libraries(radiation ${PETSC_LIBRARIES})

#set_source_files_properties(radiationPy.pyx PROPERTIES CYTHON_IS_CXX TRUE)
#cython_add_module(radiationPy radiationPy.pyx)
#target_link_libraries(radiationPy radiation physics grid geometry ${ArrayFire_LIBRARIES})
