add_library(problem ${PROBLEM_DIR}/problem.c ${PROBLEM_DIR}/problemdata.h problem.h)
target_link_libraries(problem timestepper gridzone physics geometry boundary
                      ${MATH_LIBRARIES})
