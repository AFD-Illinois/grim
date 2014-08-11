set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")
add_library(boundary boundary.c boundary.h)
target_link_libraries(boundary gridzone physics)
