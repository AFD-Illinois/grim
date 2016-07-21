add_library(problem ${PROBLEM_DIR}/problem.cpp)
add_library(params ${PROBLEM_DIR}/params.cpp)
target_link_libraries(problem params geometry)
