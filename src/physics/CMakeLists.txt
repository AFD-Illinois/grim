add_library(physics physics.cpp physics.hpp riemannsolver.cpp
            ${CMAKE_SOURCE_DIR}/problem/${PROBLEM}/params.cpp
           )
target_link_libraries(physics reconstruction)
