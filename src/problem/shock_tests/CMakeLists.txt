# Dump folder
set(OUTPUT_DIR "\"/home/mc/data/mc/PhD/grim/build/\"")
set(RESTART "OFF")
set(RESTART_FILE "\"restartfile.h5\"")
set(DUMP_FILE_PREFIX "\"data\"")
set(PROBLEM_DATA "\"${PROBLEM}/problemdata.h\"")

# Time stepping options: EXPLICIT, IMEX or IMPLICIT
set(TIME_STEPPING "EXPLICIT")
set(DT "1e-10")
set(DT_DUMP ".1")
set(START_TIME "0.")
set(FINAL_TIME "2.5")
set(COURANT "0.5")
set(MAX_DT_INCREMENT "1.3")

# Domain size. If the problem is 1D then N2 is ignored.
set(COMPUTE_DIM "1")
set(N1 "2048")
set(N2 "128")

# The entire global domain is divided into tiles that are small enough to fit
# into the cache of the compute node or an accelerator. This technique
# optimizes cache usage and makes explicit use of the deep memory hierarchies
# which are prevalent in all kinds of processors and accelerators. 
# 
# Caution: 1) The tile sizes need to divide the global domain size in each
#             direction exactly!
#          2) We want a tile size that is as large as the cache. Too large a
#             size and the code can fail silently because local memory in the
#             OpenCL specification is a different amount for different 
#             machines.
#             NEED TO PUT A CHECK FOR THIS.
set(TILE_SIZE_X1 "8")
set(TILE_SIZE_X2 "1")
set(USE_OPENMP   "NO")

# Physics variables
set(ADIABATIC_INDEX "4./3")
set(CONDUCTION "OFF")

# Geometry
set(METRIC "MINKOWSKI")
set(EPS    "1e-5")

# Domain
set(X1_A  "-5")
set(X1_B  "5.")
set(X2_A  "0.")
set(X2_B  "1.")

# Boundary conditions
set(PHYSICAL_BOUNDARY_LEFT_EDGE   "OUTFLOW")
set(PHYSICAL_BOUNDARY_RIGHT_EDGE  "OUTFLOW")
set(PHYSICAL_BOUNDARY_TOP_EDGE    "PERIODIC")
set(PHYSICAL_BOUNDARY_BOTTOM_EDGE "PERIODIC")

# Reconstrution options
# MP5, MONOTONIZED_CENTRAL or MIN_MOD
set(RECONSTRUCTION "MONOTONIZED_CENTRAL")

# Initial condition parameters 
# Komissarov shock tests for MHD
# Choose the shock test:
# FAST_SHOCK
# SLOW_SHOCK
# SWITCH_ON_SLOW
# SWITCH_OFF_FAST
# ALFVEN_WAVE
# SHOCK_TUBE_1
# SHOCK_TUBE_2
# COLLISION
set(SHOCK_TEST  "SLOW_SHOCK")

# Floor values
set(RHO_FLOOR         "1e-7")
SET(UU_FLOOR          "1e-9")
set(RHO_FLOOR_MIN     "1e-15")
set(UU_FLOOR_MIN      "1e-15")

message("")
message("#########################")
message("# Configuration options #")
message("#########################")
message("")
message("Problem          : " ${PROBLEM})
message("Output dir       : " ${OUTPUT_DIR})
message("Time stepping    : " ${TIME_STEPPING})
message("Dimensions       : " ${COMPUTE_DIM})
message("Resolution       : " ${N1})
message("Tile size        : " ${TILE_SIZE_X1})
message("Metric           : " ${METRIC})

message("")
message("##################################")
message("# Domain and boundary conditions #")
message("##################################")
message("")
message("Left edge      : " ${X1_A})
message("Right edge     : " ${X1_B})
message("Left boundary  : " ${PHYSICAL_BOUNDARY_LEFT_EDGE})
message("Right boundary : " ${PHYSICAL_BOUNDARY_RIGHT_EDGE})

message("")
message("###################")
message("# Physics options #")
message("###################")
message("")
message("Adiabatic index  : " ${ADIABATIC_INDEX})
message("Conduction       : " ${CONDUCTION})

message("")
message("##################")
message("# Reconstruction #")
message("##################")
message("")
message("Reconstruction  : " ${RECONSTRUCTION})

message("")
message("######################")
message("# Initial conditions #")
message("######################")
message("")
message("Shock test : " ${SHOCK_TEST})

message("")
message("##########")
message("# Floors #")
message("##########")
message("")
message("Density floor            : " ${RHO_FLOOR})
message("Internal energy floor    : " ${UU_FLOOR})
