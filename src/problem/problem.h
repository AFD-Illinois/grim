#ifndef GRIM_PROBLEM_H_
#define GRIM_PROBLEM_H_

#define LINEAR_MODES_1D                                   (0)
#define KOMISSAROV_FAST_SHOCK_1D                          (1)
#define KOMISSAROV_SLOW_SHOCK_1D                          (2)
#define KOMISSAROV_SWITCH_OFF_1D                          (3)      
#define KOMISSAROV_SWITCH_ON_1D                           (4)
#define KOMISSAROV_SHOCK_TUBE_1_1D                        (5)
#define KOMISSAROV_SHOCK_TUBE_2_1D                        (6)
#define KOMISSAROV_COLLISION_1D                           (7)
#define MAGNETIZED_FIELD_LOOP_ADVECTION_2D                (8)
#define ORZAG_TANG_2D                                     (9)
#define KOMISSAROV_MAGNETIZED_CYLINDRICAL_EXPLOSION_2D    (10)
#define EQUILIBRIUM_TORUS_2D                              (11)
#define MHD_TORUS_2D                                      (12)
#define MTI_TEST_2D                                       (13)
#define MHD_WITH_CONDUCTION_TORUS_2D                      (14)


#if (PROBLEM==LINEAR_MODES_1D)
#include "linear_modes_1d.h"
#elif (PROBLEM==KOMISSAROV_FAST_SHOCK_1D)
#include "komissarov_fast_shock_1d.h"
#elif (PROBLEM==KOMISSAROV_SLOW_SHOCK_1D)
#include "komissarov_slow_shock_1d.h"
#elif (PROBLEM==KOMISSAROV_SWITCH_OFF_1D)
#include "komissarov_switch_off_1d.h"
#elif (PROBLEM==KOMISSAROV_SWITCH_ON_1D)
#include "komissarov_switch_on_1d.h"
#elif (PROBLEM==KOMISSAROV_SHOCK_TUBE_1_1D)
#include "komissarov_shock_tube_1_1d.h"
#elif (PROBLEM==KOMISSAROV_SHOCK_TUBE_2_1D)
#include "komissarov_shock_tube_2_1d.h"
#elif (PROBLEM==KOMISSAROV_COLLISION_1D)
#include "komissarov_collision_1d.h"
#elif (PROBLEM==MAGNETIZED_FIELD_LOOP_ADVECTION_2D)
#include "magnetized_field_loop_advection_2d.h"
#elif (PROBLEM==ORZAG_TANG_2D)
#include "orzag_tang_2d.h"
#elif (PROBLEM==KOMISSAROV_MAGNETIZED_CYLINDRICAL_EXPLOSION_2D)
#include "komissarov_magnetized_cylindrical_explosion_2d.h"
#elif (PROBLEM==EQUILIBRIUM_TORUS_2D)
#include "equilibrium_torus_2d.h"
#elif (PROBLEM==MHD_TORUS_2D)
#include "mhd_torus_2d.h"
#elif (PROBLEM==MTI_TEST_2D)
#include "mti_test_2d.h"
#elif (PROBLEM==MHD_WITH_CONDUCTION_TORUS_2D)
#include "mhd_with_conduction_torus_2d.h"
#endif /* Include the chosen problem header */

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1]);

void applyFloor(REAL primTile[ARRAY_ARGS TILE_SIZE]);

void postStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

#endif /* GRIM_PROBLEM_H_ */
