add_library(problem ${PROBLEM_DIR}/problem.c problem.h)
target_link_libraries(problem timestepper gridzone physics geometry boundary
                      ${MATH_LIBRARIES})
