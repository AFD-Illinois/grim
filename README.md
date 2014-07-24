grim
====

General Relativistic Implicit Magnetohydrodynamics

Dependency chart:

                                                                      gridzone
                                                                         |\______________________________________________________
                                                                         |                  params                               \
                                                                         |                    |                                  |
                                                                         |                    |                                  |
                                                                         |                 geometry                       petsc  | 
                                                                          \______   ______/   |                             |    |
                                                                                 \ /          |                             |    | 
                                                                               boundary    physics                          |    |
                                                                                  |           |                             |    |
                                                                                   \_________/ \_________                   |    |
                                                                                        |               |\                 /|    |
                                                                                        |               | initialconditions |    |
                                                                  reconstruction   riemannsolver        |             _____/|    |
                                                                        \____________   |                \           /      |    |
                                                                                     \  |                 diagnostics       |    |
                                                                                    zoneresidual                            |    |
                                                                                        |                                   |    |
                                                                                        |                                   |    |
                                                                                         \              ___________________/ \___| 
                                                                                          \            /                        /
                                                                                          globalresidual                       /
                                                                                                |                             /
                                                                                                |                        globalgrid
                                                                                             timestep                   /   
                                                                                                | _____________________|
                                                                                                |/
                                                                                              grim
                                                                                          
                                                                                          
                                                                                          
                                                                                          
                                                                                          
                                                                                          
                                                                                          
