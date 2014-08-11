set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")
add_library(reconstruct reconstruct.c reconstruct.h)
target_link_libraries(reconstruct physics)
