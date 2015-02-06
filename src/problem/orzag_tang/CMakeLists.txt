# Dump folder
set(OUTPUT_DIR "\"/home/mc/data/mc/PhD/grim/build/\"")
set(RESTART "OFF")
set(RESTART_FILE "\"restartfile.h5\"")
set(DUMP_FILE_PREFIX "\"data\"")
set(PROBLEM_DATA "\"${PROBLEM}/problemdata.h\"")

# Time stepping options: EXPLICIT, IMEX or IMPLICIT
set(TIME_STEPPING "EXPLICIT")
set(DT "0.002")
set(DT_DUMP ".05")
set(START_TIME "0.")
set(FINAL_TIME "1.")
set(COURANT "0.5")
set(MAX_DT_INCREMENT "1.3")

# Domain size. If the problem is 1D then N2 is ignored.
set(COMPUTE_DIM "2")
set(N1 "280")
set(N2 "280")

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
set(TILE_SIZE_X1 "10")
set(TILE_SIZE_X2 "10")
set(USE_OPENMP   "NO")

# Physics variables
set(ADIABATIC_INDEX "5./3")
set(CONDUCTION "OFF")

# Geometry
set(METRIC "MINKOWSKI")
set(EPS    "1e-5")

# Domain
set(X1_A  "0.")
set(X1_B  "1.")
set(X2_A  "0.")
set(X2_B  "1.")

# Boundary conditions
set(PHYSICAL_BOUNDARY_LEFT_EDGE   "PERIODIC")
set(PHYSICAL_BOUNDARY_RIGHT_EDGE  "PERIODIC")
set(PHYSICAL_BOUNDARY_TOP_EDGE    "PERIODIC")
set(PHYSICAL_BOUNDARY_BOTTOM_EDGE "PERIODIC")

# Reconstrution options
# MP5, MONOTONIZED_CENTRAL or MIN_MOD
set(RECONSTRUCTION "MP5")

# Floor values
set(RHO_FLOOR "1e-5")
SET(UU_FLOOR "1e-7")

message("")
message("#########################")
message("# Configuration options #")
message("#########################")
message("")
message("Problem          : " ${PROBLEM})
message("Output dir       : " ${OUTPUT_DIR})
message("Time stepping    : " ${TIME_STEPPING})
message("Dimensions       : " ${COMPUTE_DIM})
message("Resolution       : " ${N1} " x " ${N2})
message("Tile size        : " ${TILE_SIZE_X1} " x " ${TILE_SIZE_X2})
message("Metric           : " ${METRIC})

message("")
message("##################################")
message("# Domain and boundary conditions #")
message("##################################")
message("")
message("              " ${PHYSICAL_BOUNDARY_TOP_EDGE})
message("    (" ${X1_A},${X2_B} ")+-----------+" "(" ${X1_B},${X2_B} ")" )
message(                "           |           |"       )
message(                "           |           |"       )
message("  " ${PHYSICAL_BOUNDARY_LEFT_EDGE} " |           | " ${PHYSICAL_BOUNDARY_RIGHT_EDGE})
message(                "           |           |"       )
message(                "           |           |"       )
message("    (" ${X1_A},${X2_A} ")+-----------+" "(" ${X1_B},${X2_A} ")" )
message("              " ${PHYSICAL_BOUNDARY_BOTTOM_EDGE})

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
message("##########")
message("# Floors #")
message("##########")
message("")
message("Density floor            : " ${RHO_FLOOR})
message("Internal energy floor    : " ${UU_FLOOR})
