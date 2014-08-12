set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")
add_library(riemannsolver riemannsolver.c riemannsolver.h)
target_link_libraries(riemannsolver physics geometry)
