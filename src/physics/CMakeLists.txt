add_library(physics physics.cpp physics.hpp riemannsolver.cpp
            ${CMAKE_SOURCE_DIR}/problem/${PROBLEM}/params.cpp
            ${CMAKE_SOURCE_DIR}/problem/${PROBLEM}/problem.cpp
           )
target_link_libraries(physics reconstruction ${PETSC_LIBRARIES})
