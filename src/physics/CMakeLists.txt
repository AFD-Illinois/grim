set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")
add_library(physics physics.c physics.h)
target_link_libraries(physics geometry)
