add_library(timestepper timestepper.c timestepper.h residual.c fvmfluxes.c
            diagnostics.c)
target_link_libraries(timestepper geometry gridzone boundary reconstruct
                      riemannsolver problem)
