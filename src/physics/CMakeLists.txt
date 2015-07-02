add_library(physics physics.c physics.h riemannsolver.c conduction.c viscosity.c)
target_link_libraries(physics geometry problem)
